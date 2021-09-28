
#include <numeric> // for accumulate

#include "../include/quasineutral.hxx"

Quasineutral::Quasineutral(std::string name, Options &alloptions,
                           Solver *UNUSED(solver))
    : name(name) {
  Options &options = alloptions[name];
  
  // Need to have a charge and mass
  charge = options["charge"].doc("Particle charge. electrons = -1");
  AA = options["AA"].doc("Particle atomic mass. Proton = 1");

  // Save the density
  bout::globals::dump.addRepeat(density, std::string("N") + name);
  density = 0.0;

  ASSERT0(charge != 0.0);
}

void Quasineutral::transform(Options &state) {
  // Iterate through all subsections
  Options &allspecies = state["species"];

  // Add charge density of other species
  Field3D rho = std::accumulate(
      // Iterate through species
      begin(allspecies.getChildren()), end(allspecies.getChildren()),
      // Start with no charge
      Field3D(0.0),
      [this](Field3D value,
             const std::map<std::string, Options>::value_type &name_species) {
        const Options &species = name_species.second;
        // Add other species which have density and charge
        if (name_species.first != name and species.isSet("charge") and
            species.isSet("density")) {
          // Note: Not assuming that the boundary has been set
          return value + getNoBoundary<Field3D>(species["density"]) *
                             get<BoutReal>(species["charge"]);
        }
        return value;
      });

  // Set quantites for this species
  Options &species = allspecies[name];

  // Calculate density required. Floor so that density is >= 0
  density = floor(rho / (-charge), 0.0);
  set(species["density"], density);

  set(species["charge"], charge);

  set(species["AA"], AA);
}
