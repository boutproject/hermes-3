
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

  offset = options["potential_offset"]
               .doc("Potential at which the sheath current is zero")
               .withDefault<BoutReal>(0.0);

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
}

