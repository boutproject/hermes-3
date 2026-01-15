#include <bout/fv_ops.hxx>
#include <bout/vecops.hxx>
#include <bout/constants.hxx>

#include "../include/parallel_inertia_curvature_drift.hxx"

using bout::globals::mesh;

ParallelInertiaCurvatureDrift::ParallelInertiaCurvatureDrift(std::string name, Options& alloptions,
                                   Solver* UNUSED(solver)) {

  // Get options for this component
  auto& options = alloptions[name];

  bndry_flux =
      options["bndry_flux"].doc("Allow fluxes through boundary?").withDefault<bool>(true);

  diagnose =
      options["diagnose"].doc("Output additional diagnostics?").withDefault<bool>(false);
      
  // Read curvature vector
  Curlb_B.covariant = false; // Contravariant
  if (mesh->get(Curlb_B, "bxcv")) {
    Curlb_B.x = Curlb_B.y = Curlb_B.z = 0.0;
  }

  Options& paralleltransform = Options::root()["mesh"]["paralleltransform"];
  if (paralleltransform.isSet("type") and
      paralleltransform["type"].as<std::string>() == "shifted") {
    Field2D I;
    if (mesh->get(I, "sinty")) {
      I = 0.0;
    }
    Curlb_B.z += I * Curlb_B.x;
  }

  // Normalise

  // Get the units
  const auto& units = alloptions["units"];
  BoutReal Bnorm = get<BoutReal>(units["Tesla"]);
  BoutReal Lnorm = get<BoutReal>(units["meters"]);

  Curlb_B.x /= Bnorm;
  Curlb_B.y *= SQ(Lnorm);
  Curlb_B.z *= SQ(Lnorm);

  Curlb_B *= 2. / mesh->getCoordinates()->Bxy;

  // Set drift to zero through sheath boundaries.
  // Flux through those cell faces should be set by sheath.
  for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
    Curlb_B.y(r.ind, mesh->ystart - 1) = -Curlb_B.y(r.ind, mesh->ystart);
  }
  for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    Curlb_B.y(r.ind, mesh->yend + 1) = -Curlb_B.y(r.ind, mesh->yend);
  }

  // Initialise current field
  DivJicurv = 0.0;
}

void ParallelInertiaCurvatureDrift::transform(Options& state) {
  // Iterate through all subsections
  Options& allspecies = state["species"];

  for (auto& kv : allspecies.getChildren()) {
    Options& species = allspecies[kv.first]; // Note: Need non-const

    //NOTE: Should we skip electrons?
    // if (kv.first == "e") {
    //   continue; // Skip electrons -> only ions
    // }

    if (!(species.isSet("charge") and species.isSet("velocity") and species.isSet("AA")))
      continue; // Skip, go to next species

    // Calculate parallel inertia curvature drift velocity for this species
    auto q = get<BoutReal>(species["charge"]);
    if (fabs(q) < 1e-5) {
      continue;
    }

    auto Vpar = GET_VALUE(Field3D, species["velocity"]); // Parallel velocity
    const auto AA = get<BoutReal>(species["AA"]);

    // Parallel inertia curvature drift velocity
    Vector3D vI = 0.5 * AA / q * SQ(Vpar) * Curlb_B;

    if (IS_SET(species["density"])) {
      auto N = GET_VALUE(Field3D, species["density"]);

      Field3D density_source = FV::Div_f_v(N, vI, bndry_flux);

      subtract(species["density_source"], density_source);


      Vector3D Jicurv = N * vI * q; // Parallel inertia curvature current for this species

      // Divergence of current in vorticity equation
      DivJicurv += Div(Jicurv);

    }

    if (IS_SET(species["pressure"])) {
      auto P = get<Field3D>(species["pressure"]);

      Field3D energy_source = (5. / 2) * FV::Div_f_v(P, vI, bndry_flux);
      subtract(species["energy_source"],  energy_source );
    }

    if (IS_SET(species["momentum"])) {
      auto NV = get<Field3D>(species["momentum"]);
      Field3D momentum_source = FV::Div_f_v(NV, vI, bndry_flux);
      subtract(species["momentum_source"], momentum_source);
    }

  }

  add(state["fields"]["DivJextra"], DivJicurv);

}

void ParallelInertiaCurvatureDrift::outputVars(Options& state) {
  AUTO_TRACE();

  if (diagnose) {

    // Normalisations
    auto Nnorm = get<BoutReal>(state["Nnorm"]);
    auto Omega_ci = get<BoutReal>(state["Omega_ci"]);

    set_with_attrs(state["DivJicurv"], DivJicurv,
                    {{"time_dimension", "t"},
                    {"units", "A m^-3"},
                    {"conversion", SI::qe * Nnorm * Omega_ci},
                    {"long_name", "Divergence of parallel inertia curvature current"},
                    {"source", "parallel_inertia_curvature_drift"}});

  }
}