
#include <bout/constants.hxx>

#include <bout/mesh.hxx> 

#include "../include/isothermal.hxx"

using bout::globals::mesh;

Isothermal::Isothermal(std::string name, Options &alloptions,
                       Solver *UNUSED(solver))
    : name(name) {
  AUTO_TRACE();
  Options& options = alloptions[name];

  auto Tnorm = get<BoutReal>(alloptions["units"]["eV"]);
  T = options["temperature"].doc("Constant temperature [eV]").as<BoutReal>()
      / Tnorm; // Normalise

  diagnose = options["diagnose"]
    .doc("Save additional output diagnostics")
    .withDefault<bool>(false);

}

void Isothermal::transform(Options &state) {
  AUTO_TRACE();

  Options& species = state["species"][name];

  T_f3d = T;
  if (T_f3d.isFci()) {
    T_f3d.applyBoundary("neumann_o2");
    mesh->communicate(T_f3d);
    T_f3d.applyParallelBoundary("parallel_neumann_o2");
  }
  
  set(species["temperature"], T_f3d);

  // If density is set, also set pressure
  if (isSetFinalNoBoundary(species["density"])) {
    // Note: The boundary of N may not be set yet
    auto N = GET_NOBOUNDARY(Field3D, species["density"]);
    P = N * T_f3d;
    set(species["pressure"], P);
  }
}

void Isothermal::outputVars(Options& state) {
  AUTO_TRACE();
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Nnorm = get<BoutReal>(state["Nnorm"]);

  // Save the temperature to the output files
  set_with_attrs(state[std::string("T") + name], T,
                 {{"units", "eV"},
                  {"conversion", Tnorm},
                  {"long_name", name + " temperature"},
                  {"standard_name", "temperature"},
                  {"species", name},
                  {"source", "isothermal"}});

  if (diagnose) {
    // Save pressure as time-varying field
    set_with_attrs(state[std::string("P") + name], P,
                   {{"time_dimension", "t"},
                    {"units", "Pa"},
                    {"conversion", SI::qe * Tnorm * Nnorm},
                    {"long_name", name + " pressure"},
                    {"standard_name", "pressure"},
                    {"species", name},
                    {"source", "isothermal"}});
   }
}
