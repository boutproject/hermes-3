
#include "../include/sheath_closure.hxx"

SheathClosure::SheathClosure(std::string name, Options &alloptions, Solver *) {
  Options& options = alloptions[name];

  BoutReal Lnorm = alloptions["units"]["meters"]; // Length normalisation factor
  
  L_par = options["connection_length"].as<BoutReal>() / Lnorm;

  output.write("L_par = %e\n", L_par);
}

void SheathClosure::transform(Options &state) {
  AUTO_TRACE();
  
  // Get electrostatic potential
  auto phi = get<Field3D>(state["fields"]["phi"]);

  auto n = get<Field3D>(state["species"]["e"]["density"]);
  
  add(state["fields"]["DivJextra"], // Used in vorticity
      n * phi / L_par);

  add(state["species"]["e"]["density_source"],
      n * phi / L_par);
}

