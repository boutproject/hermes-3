#include "../include/external_apar.hxx"
#include "../include/component.hxx"
#include "../include/guarded_options.hxx"
#include "../include/permissions.hxx"

#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/globals.hxx>
#include <bout/mesh.hxx>
#include <bout/options.hxx>

#include <string>

ExternalApar::ExternalApar(std::string name, Options& alloptions,
                           [[maybe_unused]] Solver* solver)
    : NamedComponent(name, {readWrite("fields:Apar_flutter")}) {

  auto& options = alloptions[name];
  const std::string apar_name = options["apar_name"]
                                    .doc("Name of the Apar field in the mesh file")
                                    .withDefault("external_apar");

  // Read a 3D field from the input e.g. mesh file
  // Store in member variable
  if (bout::globals::mesh->get(external_apar, apar_name) != 0) {
    throw BoutException(
        "Couldn't read `{}` in external_apar. See the `apar_name` option.", apar_name);
  }

  // Normalise and scale
  const Options& units = alloptions["units"];
  const BoutReal rho_s0 = units["meters"];
  const BoutReal Bnorm = units["Tesla"];

  const BoutReal scale =
      options["scale"].doc("Multiply external Apar by this factor").withDefault(1.0);

  external_apar *= scale / (Bnorm * rho_s0);
}

void ExternalApar::transform_impl(GuardedOptions& state) {
  // Add the field member variable to Apar_flutter
  add(state["fields"]["Apar_flutter"], external_apar);
}

void ExternalApar::outputVars(Options& state) {
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
