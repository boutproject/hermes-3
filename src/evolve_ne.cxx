
#include <bout/fv_ops.hxx>
#include <derivs.hxx>
#include <difops.hxx>
#include <bout/constants.hxx>

#include "../include/evolve_ne.hxx"
#include "../include/div_ops.hxx"

using bout::globals::mesh;

EvolveNe::EvolveNe(std::string name, Options &alloptions, Solver *solver) {
  AUTO_TRACE();
  
  // Evolve the electron density in time
  solver->add(Ne, "Ne");

  auto& options = alloptions[name];

  ne_bndry_flux = options["bndry_flux"]
                      .doc("Allow flows through radial boundaries")
                      .withDefault<bool>(true);

  poloidal_flows = options["poloidal_flows"]
                       .doc("Include poloidal ExB flow")
                       .withDefault<bool>(true);

  anomalous_D = options["anomalous_D"]
                    .doc("Anomalous diffusion (axisymmetric). <0 => off")
                    .withDefault(-1.0);

  low_n_diffuse = options["low_n_diffuse"]
                      .doc("Parallel diffusion at low density")
                      .withDefault<bool>(false);

  low_n_diffuse_perp = options["low_n_diffuse_perp"]
                           .doc("Perpendicular diffusion at low density")
                           .withDefault<bool>(false);
  
  hyper_z = options["hyper_z"].doc("Hyper-diffusion in Z").withDefault<bool>(false);
}

void EvolveNe::transform(Options &state) {
  AUTO_TRACE();
  mesh->communicate(Ne);
  
  auto& electrons = state["species"]["e"];
  set(electrons["density"], Ne);
  set(electrons["AA"], SI::Me / SI::Mp); // Atomic mass
  set(electrons["charge"], -1.0);
}

void EvolveNe::finally(const Options &state) {
  AUTO_TRACE();

  // Get the coordinate system
  auto coord = Ne.getCoordinates();
  
  auto& electrons = state["species"]["e"];

  if (state.isSection("fields") and state["fields"].isSet("phi")) {
    // Electrostatic potential set -> include ExB flow

    Field3D phi = get<Field3D>(state["fields"]["phi"]);

    ddt(Ne) = -Div_n_bxGrad_f_B_XPPM(Ne, phi, ne_bndry_flux, poloidal_flows,
                                     true); // ExB drift
  } else {
    ddt(Ne) = 0.0;
  }

  if (electrons.isSet("velocity")) {
    // Parallel velocity set
    Field3D Ve = get<Field3D>(electrons["velocity"]);
    
    // Typical wave speed used for numerical diffusion
    Field3D sound_speed;
    if (state.isSet("sound_speed")) {
      Field3D sound_speed = get<Field3D>(state["sound_speed"]);
    } else {
      Field3D Te = get<Field3D>(electrons["temperature"]);
      sound_speed = sqrt(Te);
    }

    if (state.isSection("fields") and state["fields"].isSet("phi")) {
      // Parallel wave speed increased to electron sound speed
      // since electrostatic & electromagnetic waves are supported
      ddt(Ne) -= FV::Div_par(Ne, Ve, sqrt(SI::Me / SI::Mp) * sound_speed);
    } else {
      // Parallel wave speed is ion sound speed
      ddt(Ne) -= FV::Div_par(Ne, Ve, sound_speed);
    }
  }
  
  if (anomalous_D > 0.0) {
    ddt(Ne) += FV::Div_a_Laplace_perp(anomalous_D, DC(Ne));
  }

  if (low_n_diffuse) {
    // Diffusion which kicks in at very low density, in order to
    // help prevent negative density regions
    ddt(Ne) += FV::Div_par_K_Grad_par(SQ(coord->dy) * coord->g_22 * 1e-4 / Ne, Ne);
  }
  if (low_n_diffuse_perp) {
    ddt(Ne) += Div_Perp_Lap_FV_Index(1e-4 / Ne, Ne, ne_bndry_flux);
  }

  if (hyper_z > 0.) {
    ddt(Ne) -= hyper_z * SQ(SQ(coord->dz)) * D4DZ4(Ne);
  }

  if (electrons.isSet("density_source")) {
    ddt(Ne) += get<Field3D>(electrons["density_source"]);
  }
}
