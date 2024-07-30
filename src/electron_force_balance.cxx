
#include <bout/difops.hxx>

#include "../include/electron_force_balance.hxx"

using bout::globals::mesh;

void ElectronForceBalance::transform(Options &state) {
  AUTO_TRACE();

  if (IS_SET(state["fields"]["phi"])) {
    // Here we use electron force balance to calculate the parallel electric field
    // rather than the electrostatic potential
    throw BoutException("Cannot calculate potential and use electron force balance\n");
  }

  // Get the electron pressure
  // Note: The pressure boundary can be set in sheath boundary condition
  //       which depends on the electron velocity being set here first.
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
  const Field3D Epar = force_density / floor(Ne, 1e-5);

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
