
#include "../include/isothermal.hxx"
#include <bout/constants.hxx>
#include <bout/mesh.hxx>
#include <bout/globalfield.hxx>
using bout::globals::mesh;

Isothermal::Isothermal(std::string name, Options& alloptions, Solver* UNUSED(solver))
    : NamedComponent(
          name, {readIfSet(fmt::format("species:{}:density", name), Regions::Interior),
                 readWrite(fmt::format("species:{}:temperature", name)),
                 // FIXME: This is only written if density is set
                 readWrite(fmt::format("species:{}:pressure", name))}) {
  Options& options = alloptions[name];

  auto Tnorm = get<BoutReal>(alloptions["units"]["eV"]);
  T = options["temperature"].doc("Constant temperature [eV]").as<BoutReal>()
      / Tnorm; // Normalise
  temperature_3D = options["temperature_3D"]
    .doc("Should the temperature be a 3D field (true) or a constant (false)?")
    .withDefault<bool>(true);
  diagnose = options["diagnose"]
    .doc("Save additional output diagnostics")
    .withDefault<bool>(false);
  dipole_scaling = options["dipole_scaling"]
                 .doc("Apply dipole scaling ~ B to the temperature")
                 .withDefault<bool>(false);

  if (dipole_scaling && !temperature_3D) {
    throw std::runtime_error("Isothermal: dipole_scaling requires temperature_3D = true");
  }

  T2D.allocate();
  B2D.allocate();
  B2D_edge.allocate();
  if (dipole_scaling) {
    
     // Dipole scaling is only implemented for 3D temperature profile only at the moment
    B2D_edge = 0;
  mesh->communicate(B2D_edge);
    mesh->get(B2D_edge, "B_edge", 0.0, true);
    T2D.allocate();
    mesh->communicate(B2D_edge);
      for (int i = mesh->xstart-2; i <= mesh->xend+2; i++) {
        for (int j = mesh->ystart; j <= mesh->yend; j++) {
          T2D(i, j) = T * mesh->getCoordinates()->Bxy(i, j) / B2D_edge(i, j);
        }
      }
      mesh->communicate(T2D);
  }

}


void Isothermal::transform_impl(GuardedOptions& state) {

  GuardedOptions species = state["species"][objectName()];
  Field3D T_= T;
  if (temperature_3D)
  {
  
  if (dipole_scaling) {
    // Set the temperature to a dipole profile, using the pressure if it's set, otherwise using a fixed temperature
    set(species["temperature"], T2D);T_ = T2D;
  }
  else{T_ = T;}
  }

  
  set(species["temperature"], T_);

  // If density is set, also set pressure
  if (isSetFinalNoBoundary(species["density"])) {
    // Note: The boundary of N may not be set yet
    auto N = GET_NOBOUNDARY(Field3D, species["density"]);
      P = N * T_;
    set(species["pressure"], P);
  }
}

void Isothermal::outputVars(Options& state) {
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  const std::string& name = objectName();

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
                    if (dipole_scaling) {
                      set_with_attrs(state[std::string("T2D") + name], T2D,
                                     {{"time_dimension", "t"},
                                      {"units", "eV"},
                                      {"conversion", Tnorm},
                                      {"long_name", name + " temperature with dipole scaling"},
                                      {"standard_name", "temperature with dipole scaling"},
                                      {"species", name},
                                      {"source", "isothermal"}});
                    
    set_with_attrs(state[std::string("B_edge")], B2D_edge,
                   {{"units", "T"}, {"long_name", "Magnetic field at edge"}, {"source", "isothermal"}});}
   }
}
