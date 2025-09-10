#include "../include/external_apar.hxx"

#include <bout/constants.hxx>
#include <bout/mesh.hxx>

ExternalApar::ExternalApar(std::string name, Options& alloptions,
                           Solver* UNUSED(solver)) {
  AUTO_TRACE();

  // Read a 3D field from the input e.g. mesh file
  // Store in member variable
  bout::globals::mesh->get(external_apar, "external_apar");

  // Normalise
  const Options& units = alloptions["units"];
  const BoutReal rho_s0 = units["meters"];
  const BoutReal Bnorm = units["Tesla"];
  external_apar /= Bnorm * rho_s0;
}

void ExternalApar::transform(Options& state) {
  AUTO_TRACE();

  // Add the field member variable to Apar_flutter
  add(state["fields"]["Apar_flutter"], external_apar);
}

void ExternalApar::outputVars(Options& state) {
  AUTO_TRACE();

  // Normalisations
  auto Bnorm = get<BoutReal>(state["Bnorm"]);
  auto rho_s0 = get<BoutReal>(state["rho_s0"]);

  set_with_attrs(state["external_apar"], external_apar,
                 {{"units", "Tm"},
                  {"conversion", Bnorm * rho_s0},
                  {"standard_name", "External Apar"},
                  {"long_name", "External Apar"},
                  {"source", "external_apar"}});
}
