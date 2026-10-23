#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/globals.hxx>
#include <bout/output.hxx>
#include <bout/solver.hxx>

#include "../include/hermes_utils.hxx"
#include "../include/neutral_radial.hxx"

using bout::globals::mesh;

NeutralRadial::NeutralRadial(const std::string& name, Options& alloptions, Solver* solver)
    : NamedComponent(name, {readWrite("species:{name}:{outputs}")}), name(name) {
  const Options& units = alloptions["units"];
  const BoutReal meters = units["meters"];
  const BoutReal seconds = units["seconds"];
  const BoutReal Nnorm = units["inv_meters_cubed"];
  const BoutReal Tnorm = units["eV"];

  // Radial derivatives need x guard cells
  ASSERT0(mesh->xstart > 0);

  auto& options = alloptions[name];

  solver->add(Nn, std::string("N") + name);

  AA = options["AA"].doc("Particle atomic mass. Proton = 1").withDefault(1.0);
  Tn = options["temperature"].doc("Fixed neutral temperature [eV]").withDefault(3.0)
       / Tnorm;
  density_floor = options["density_floor"]
                      .doc("Minimum density, normalised units")
                      .withDefault(1e-8);
  neutral_lmax = options["neutral_lmax"]
                     .doc("Maximum neutral mean free path [m], used as a "
                          "collisionality floor")
                     .withDefault(1.0)
                 / meters;
  collisionality_override =
      options["collisionality_override"]
          .doc("Fixed neutral collision frequency in Dnn (normalised). < 0: off")
          .withDefault(-1.0);
  diffusion_limit = options["diffusion_limit"]
                        .doc("Upper limit on the diffusion coefficient [m^2/s]. < 0: off")
                        .withDefault(-1.0)
                    / (meters * meters / seconds);
  diffusion_collisions_mode = options["diffusion_collisions_mode"]
                                  .doc("afn: CX, IZ and NN collisions; multispecies: "
                                       "all collisions and CX")
                                  .withDefault<std::string>("afn");
  diagnose =
      options["diagnose"].doc("Save additional diagnostics?").withDefault<bool>(false);

  // Density source, from the mesh and/or the input [m^-3 s^-1]
  density_source = 0.0;
  mesh->get(density_source, std::string("N") + name + "_src");
  density_source = alloptions[std::string("N") + name]["source"]
                       .doc("Source term in ddt(N" + name + "). Units [m^-3/s]")
                       .withDefault(density_source)
                   / (Nnorm / seconds);

  // Dirichlet (zero) diffusivity at the walls: no diffusive flux through x boundaries
  alloptions[std::string("Dnn") + name]["bndry_all"] =
      alloptions[std::string("Dnn") + name]["bndry_all"].withDefault("dirichlet");
  alloptions[std::string("N") + name]["bndry_all"] =
      alloptions[std::string("N") + name]["bndry_all"].withDefault("neumann");
  Dnn.setBoundary(std::string("Dnn") + name);
  Nn.setBoundary(std::string("N") + name);

  substitutePermissions("name", {name});
  substitutePermissions(
      "outputs", {"AA", "density", "pressure", "temperature", "momentum", "velocity"});
}

void NeutralRadial::transform_impl(GuardedOptions& state) {
  mesh->communicate(Nn);
  Nn = floor(Nn, 0.0);

  // No parallel flow: velocity and momentum are zero, so reactions see static neutrals
  auto localstate = state["species"][name];
  set(localstate["density"], Nn);
  set(localstate["AA"], AA);
  set(localstate["temperature"], Field3D{Tn});
  set(localstate["pressure"], Field3D{Nn * Tn});
  set(localstate["velocity"], Field3D{0.0});
  set(localstate["momentum"], Field3D{0.0});
}

void NeutralRadial::finally(const Options& state) {
  const auto& localstate = state["species"][name];
  Nn = get<Field3D>(localstate["density"]);

  // Diffusion coefficient Dnn = (Tn/AA) / nu, as in neutral_mixed, with the collision
  // frequency floored by the maximum mean free path
  const BoutReal nu_pseudo_mfp = sqrt(Tn / AA) / neutral_lmax;
  Field3D Dnn_unlimited;
  if (collisionality_override > 0.0) {
    Dnn_unlimited = (Tn / AA) / collisionality_override;
  } else if (localstate.isSet("collision_frequency")) {
    if (collision_names.empty()) { // Select the collisions once
      for (const auto& collision : localstate["collision_frequencies"].getChildren()) {
        const std::string collision_name = collision.second.name();
        const bool use =
            (diffusion_collisions_mode == "afn")
                ? (collisionSpeciesMatch(collision_name, name, "+", "cx", "partial")
                   or collisionSpeciesMatch(collision_name, name, "+", "iz", "partial")
                   or collisionSpeciesMatch(collision_name, name, name, "coll", "exact"))
                : (collisionSpeciesMatch(collision_name, name, "", "cx", "partial")
                   or collisionSpeciesMatch(collision_name, name, "", "coll", "partial"));
        if (use) {
          collision_names.push_back(collision_name);
        }
      }
      if (diffusion_collisions_mode != "afn"
          and diffusion_collisions_mode != "multispecies") {
        throw BoutException("diffusion_collisions_mode for {} must be afn or "
                            "multispecies",
                            name);
      }
      if (collision_names.empty()) {
        throw BoutException("No collisions found for {} in neutral_radial", name);
      }
    }
    Field3D nu{zeroFrom(Nn)};
    for (const auto& collision_name : collision_names) {
      nu += get<Field3D>(localstate["collision_frequencies"][collision_name]);
    }
    // Smooth floor at the mean-free-path collisionality
    const Field3D nu_floored = nu + nu_pseudo_mfp * exp(-nu / nu_pseudo_mfp);
    Dnn_unlimited = (Tn / AA) / nu_floored;
  } else {
    Dnn_unlimited = Field3D{(Tn / AA) / nu_pseudo_mfp};
  }

  Dnn = Dnn_unlimited;
  if (diffusion_limit > 0.0) {
    // Harmonic blend with the upper limit, as neutral_mixed
    Dnn = Dnn_unlimited * diffusion_limit / (Dnn_unlimited + diffusion_limit);
  }
  mesh->communicate(Dnn);
  Dnn.applyBoundary();

  // Radial diffusion only: x fluxes of Div_a_Grad_perp, no y fluxes
  const Coordinates* coord = mesh->getCoordinates();
  Field3D div_radial{zeroFrom(Nn)};
  flow_xlow = zeroFrom(Nn);
  for (int i = mesh->xstart - 1; i <= mesh->xend; i++) {
    for (int j = mesh->ystart; j <= mesh->yend; j++) {
      for (int k = 0; k < mesh->LocalNz; k++) {
        const BoutReal fout =
            0.5 * (Dnn(i, j, k) + Dnn(i + 1, j, k))
            * (coord->J(i, j, k) * coord->g11(i, j, k)
               + coord->J(i + 1, j, k) * coord->g11(i + 1, j, k))
            * (Nn(i + 1, j, k) - Nn(i, j, k))
            / (coord->dx(i, j, k) + coord->dx(i + 1, j, k));
        div_radial(i, j, k) += fout / (coord->dx(i, j, k) * coord->J(i, j, k));
        div_radial(i + 1, j, k) -=
            fout / (coord->dx(i + 1, j, k) * coord->J(i + 1, j, k));
        // Flow is positive in the positive x direction
        flow_xlow(i + 1, j, k) = -fout * coord->dy(i, j, k) * coord->dz(i, j, k);
      }
    }
  }

  Sn = density_source;
  if (localstate.isSet("density_source")) {
    Sn += get<Field3D>(localstate["density_source"]);
  }
  ddt(Nn) = div_radial + Sn;
}

void NeutralRadial::outputVars(Options& state) {
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto Cs0 = get<BoutReal>(state["Cs0"]);
  auto rho_s0 = get<BoutReal>(state["rho_s0"]);

  state[std::string("N") + name].setAttributes({{"time_dimension", "t"},
                                                {"units", "m^-3"},
                                                {"conversion", Nnorm},
                                                {"standard_name", "density"},
                                                {"long_name", name + " number density"},
                                                {"species", name},
                                                {"source", "neutral_radial"}});
  if (diagnose) {
    set_with_attrs(state[std::string("Dnn") + name], Dnn,
                   {{"time_dimension", "t"},
                    {"units", "m^2/s"},
                    {"conversion", Cs0 * Cs0 / Omega_ci},
                    {"standard_name", "diffusion coefficient"},
                    {"long_name", name + " radial diffusion coefficient"},
                    {"source", "neutral_radial"}});
    set_with_attrs(state[std::string("S") + name + "_src"], Sn,
                   {{"time_dimension", "t"},
                    {"units", "m^-3 s^-1"},
                    {"conversion", Nnorm * Omega_ci},
                    {"standard_name", "density source"},
                    {"long_name", name + " number density source"},
                    {"species", name},
                    {"source", "neutral_radial"}});
    set_with_attrs(state[std::string("pf") + name + "_radial_xlow"], flow_xlow,
                   {{"time_dimension", "t"},
                    {"units", "s^-1"},
                    {"conversion", rho_s0 * SQ(rho_s0) * Nnorm * Omega_ci},
                    {"standard_name", "particle flow"},
                    {"long_name", name + " radial particle flow in X"},
                    {"species", name},
                    {"source", "neutral_radial"}});
  }
}
