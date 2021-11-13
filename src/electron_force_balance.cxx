
#include <difops.hxx>

#include "../include/electron_force_balance.hxx"

void ElectronForceBalance::transform(Options &state) {
  AUTO_TRACE();

  // Get the electron pressure
  // Note: The pressure boundary can be set in sheath boundary condition
  //       which depends on the electron velocity being set here first.
  Options& electrons = state["species"]["e"];
  Field3D Pe = getNoBoundary<Field3D>(electrons["pressure"]);
  Pe.applyBoundary("neumann");

  Field3D Ne = getNoBoundary<Field3D>(electrons["density"]);

  ASSERT1(get<BoutReal>(electrons["charge"]) == -1.0);

  // Force balance, E = (-∇p + F) / n
  Field3D force_density = - Grad_par(Pe);

  if (isSetFinal(electrons["momentum_source"], "zero_current")) {
    // Balance other forces from e.g. collisions
    // Note: marked as final so can't be set later
    force_density += get<Field3D>(electrons["momentum_source"]);
  }
  const Field3D Epar = force_density / Ne;

  // Now calculate forces on other species
  Options& allspecies = state["species"];
  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons
    }
    Options& species = allspecies[kv.first]; // Note: Need non-const

    if (!(species.isSet("density") and species.isSet("charge"))) {
      continue; // Needs both density and charge to experience a force
    }

    const Field3D N = getNoBoundary<Field3D>(species["density"]);
    const BoutReal charge = get<BoutReal>(species["charge"]);

    add(species["momentum_source"],
        charge * N * Epar);
  }
}
