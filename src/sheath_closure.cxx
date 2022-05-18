
#include "../include/sheath_closure.hxx"

SheathClosure::SheathClosure(std::string name, Options &alloptions, Solver *) {
  Options& options = alloptions[name];

  BoutReal Lnorm = alloptions["units"]["meters"]; // Length normalisation factor

  L_par = options["connection_length"]
              .doc("Field-line connection length in meters")
              .as<BoutReal>() /
          Lnorm;

  sheath_gamma = options["sheath_gamma"]
          .doc("Sheath heat transmission coefficient (dimensionless)")
          .withDefault<BoutReal>(6.5);

  sheath_gamma_ions = options["sheath_gamma_ions"]
          .doc("Sheath heat transmission coefficient for ions (dimensionless)")
          .withDefault<BoutReal>(2.0); // Value suggested by Stangeby's book, between eqs.
                                       // (2.92) and (2.93)

  offset = options["potential_offset"]
               .doc("Potential at which the sheath current is zero")
               .withDefault<BoutReal>(0.0);

  sinks = options["sinks"]
               .doc("Include sinks of density and energy?")
               .withDefault<bool>(false);

  output.write("\tL_par = {:e} (normalised)\n", L_par);
}

void SheathClosure::transform(Options &state) {
  AUTO_TRACE();
  
  // Get electrostatic potential
  auto phi = get<Field3D>(state["fields"]["phi"]);

  auto& electrons = state["species"]["e"];
  
  // Electron density
  auto n = get<Field3D>(electrons["density"]);

  // Divergence of current through the sheath
  Field3D DivJsh = n * (phi - offset) / L_par;
  
  add(state["fields"]["DivJextra"], // Used in vorticity
      DivJsh);

  add(electrons["density_source"], DivJsh);

  // Electron heat conduction
  if (electrons.isSet("temperature")) {
    // Assume attached, sheath-limited regime
    // Sheath heat transmission gamma * n * T * cs

    auto Te = get<Field3D>(electrons["temperature"]);
    
    Field3D qsheath = floor(sheath_gamma * n * Te * sqrt(Te), 0.0);
    
    subtract(electrons["energy_source"], qsheath / L_par);
  }

  if (sinks) {
    // Add sinks of density and energy. Allows source-driven 2d turbulence simulations.

    // Calculate a sound speed from total pressure and total mass density.
    // [Not sure if this is correct or the best thing in general, but reduces to the
    // standard Bohm boundary conditions for a pure, hydrogenic plasma.]
    Field3D P_total = 0.0;
    Field3D rho_total = 0.0; // mass density
    Options& allspecies = state["species"];
    for (auto& kv : allspecies.getChildren()) {
      Options& species = allspecies[kv.first];

      const BoutReal A = get<BoutReal>(species["AA"]);
      Field3D Ns = get<Field3D>(species["density"]);
      Field3D Ts = get<Field3D>(species["temperature"]);

      P_total += Ns * Ts;
      rho_total += A * Ns;
    }

    Field3D c_s = sqrt(P_total / rho_total);

    for (auto& kv : allspecies.getChildren()) {
      Options& species = allspecies[kv.first];
      Field3D Ns = get<Field3D>(species["density"]);

      Field3D sheath_flux = floor(Ns * c_s, 0.0);
      subtract(species["density_source"], sheath_flux / L_par);

      if (kv.first != "e") {
        // Electron sheath heat flux has different gamma-factor than ions, and was already
        // handled above

        auto Ts = get<Field3D>(species["temperature"]);

        Field3D qsheath = floor(sheath_gamma_ions * Ts * sheath_flux, 0.0);

        subtract(species["energy_source"], qsheath / L_par);
      }
    }
  }
}
