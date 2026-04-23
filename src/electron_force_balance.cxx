
#include <bout/difops.hxx>

#include "../include/electron_force_balance.hxx"
#include <bout/immersed_boundary.hxx>

using bout::globals::mesh;

void ElectronForceBalance::transform(Options &state) {
  AUTO_TRACE();

  if (IS_SET(state["fields"]["phi"])) {
    // Here we use electron force balance to calculate the parallel electric field
    // rather than the electrostatic potential
    throw BoutException("Cannot calculate potential and use electron force balance\n");
  }

  // Get the electron pressure, with boundary condition applied
  Options& electrons = state["species"]["e"];
  Field3D Pe = GET_VALUE(Field3D, electrons["pressure"]);
  Field3D Ne = GET_NOBOUNDARY(Field3D, electrons["density"]);

  ASSERT1(get<BoutReal>(electrons["charge"]) == -1.0);

  // Force balance, E = (-∇p + F) / n
  Field3D force_density = - Grad_par(Pe);

  if (IS_SET(electrons["momentum_source"])) {
    // Balance other forces from e.g. collisions
    // Note: marked as final so can't be changed later
    force_density += GET_VALUE(Field3D, electrons["momentum_source"]);
  }
  // Parallel electric field. Stored as a private member variable
  // and may be saved as an output diagnostic
  Field3D Nfloor = copy(Ne);
  if (immBndry) {
    //IB_TODO: Complex logic for / operator and late Epar setup.
    //IB_TODO: Just floor, don't need to SetBoundary again because Epar will below.
    immBndry->FloorField(Nfloor, 1e-5);
    Epar = force_density;
    BOUT_FOR(i, mesh->getRegion3D("RGN_NO_IMM_BNDRY")) {
      Epar[i] /= Nfloor[i];
    }
    Epar.name = "Epar";
    if (!immBndry->CheckFieldSetUp(Epar.name)) {
      immBndry->FieldSetup(Epar);
    }
    immBndry->SetBoundary(Epar);
  } else {
    Nfloor = floor(Ne, 1e-5);
    Epar = force_density / Nfloor;
  }

  // Now calculate forces on other species
  Options& allspecies = state["species"];
  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons
    }
    Options& species = allspecies[kv.first]; // Note: Need non-const

    if (!(IS_SET(species["density"]) and IS_SET(species["charge"]))) {
      continue; // Needs both density and charge to experience a force
    }

    const Field3D N = GET_NOBOUNDARY(Field3D, species["density"]);
    const BoutReal charge = get<BoutReal>(species["charge"]);

    add(species["momentum_source"],
        charge * N * Epar);
  }
}

void ElectronForceBalance::outputVars(Options& state) {
  if (!diagnose) {
    return;
  }
  AUTO_TRACE();
  // Get normalisations
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Lnorm = get<BoutReal>(state["rho_s0"]);

  set_with_attrs(state["Epar"], Epar,
                 {{"time_dimension", "t"},
                  {"units", "V / m"},
                  {"conversion", Tnorm / Lnorm},
                  {"standard_name", "Epar"},
                  {"long_name", "Parallel electric field"},
                  {"source", "electron_force_balance"}});
}
