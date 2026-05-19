
#include <bout/constants.hxx>

#include <bout/mesh.hxx> 

#include "../include/adhoc_potential.hxx"

using bout::globals::mesh;

AdhocPotential::AdhocPotential(std::string name, Options &alloptions, Solver *UNUSED(solver)){
  AUTO_TRACE();
  Options& options = alloptions[name];
  lambda =
      options["lambda"].doc("Sets the coefficient for phi = lambda * Te, default is 3.0, if negative calculate it urself").withDefault<BoutReal>(3.0);
  
}

void AdhocPotential::transform(Options &state) {
  AUTO_TRACE();
  auto& fields = state["fields"];
  Options& allspecies = state["species"];

  Options& electrons = allspecies["e"];

  if (!electrons.isSet("temperature")) {
    throw BoutException("Adhoc potential should be calculated, but can not find electron temperature");
  }
  
  //T_e = toFieldAligned(floor(GET_NOBOUNDARY(Field3D, electrons["temperature"]), 0.0));
  T_e = getNoBoundary<Field3D>(electrons["temperature"]);
  if (lambda > 0.0) {
    phi_adhoc = lambda * T_e;
    phi_adhoc.applyBoundary("neumann");
    mesh->communicate(phi_adhoc);
    phi_adhoc.applyParallelBoundary("parallel_neumann_o1");
  }

  set(fields["phi"], phi_adhoc);
  

}

void AdhocPotential::outputVars(Options& state) {
  AUTO_TRACE();
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Lnorm = get<BoutReal>(state["rho_s0"]);
    set_with_attrs(state["phi"], phi_adhoc,
                 {{"time_dimension", "t"},
                  {"units", "V"},
                  {"conversion", Tnorm},
                  {"standard_name", "potential"},
                  {"long_name", "plasma potential"},
                  {"source", "adhoc_potential"}});
  
}
