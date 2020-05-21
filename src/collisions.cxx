#include <iterator>

#include <bout/constants.hxx>

#include "../include/collisions.hxx"

Collisions::Collisions(std::string name, Options& alloptions, Solver*) {
  const Options& units = alloptions["units"];

  // Normalisations
  Tnorm = units["eV"];
  Nnorm = units["inv_meters_cubed"];
  rho_s0 = units["meters"];
}

void Collisions::transform(Options &state) {
  
  Options& allspecies = state["species"];

  // Treat electron collisions specially
  // electron-ion and electron-neutral collisions
  
  Options& electrons = allspecies["e"];
  const Field3D Te = get<Field3D>(electrons["temperature"]);
  
  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      // electron-electron collisions
      
      continue;
    }
    
    Options& species = allspecies[kv.first]; // Note: Need non-const

    if (species.isSet("charge")) {
      // electron-charged ion collisions
      
      
    } else {
      ////////////////////////////////////
      // electron-neutral collisions

      // Neutral density
      Field3D Nn = get<Field3D>(species["density"]);
      // Atomic mass
      BoutReal AA = get<BoutReal>(species["AA"]);
      
      BoutReal a0 = PI*SQ(5.29e-11); // Cross-section [m^2]
      
      // Electron thermal speed
      Field3D vth_e = sqrt((SI::Mp / SI::Me) * Te);
      
      // Electron-neutral collision rate
      Field3D nu_en = vth_e * Nnorm * Nn * a0 * rho_s0;

      add(electrons["collision_rate"], nu_en);

      // Elastic collisions differ by mass ratio
      add(species["collision_rate"], nu_en * (SI::Me / (SI::Mp * AA)));
    }
  }

  // Iterate through other species
  // To avoid double counting, this needs to iterate over pairs
  // i.e. the diagonal and above
  //
  // Iterators kv1 and kv2 over the species map
  // 
  //               kv2 -> 
  //             species1  species2  species3
  // kv1   species1     X         X         X
  //  ||   species2               X         X
  //  \/   species3                         X
  //
  const std::map<std::string, Options>& children = allspecies.getChildren();
  for (auto kv1 = std::begin(children); kv1 != std::end(children); ++kv1) {
    if (kv1->first == "e")
      continue; // Skip electrons
    
    Options& species1 = allspecies[kv1->first];

    // If temperature isn't set, assume zero
    const Field3D temperature1 =
        species1.isSet("temperature") ? get<Field3D>(species1["temperature"]) : 0.0;
    const BoutReal mass1 = get<BoutReal>(species1["AA"]) * SI::Mp; // in Kg
      
    if (species1.isSet("charge")) {
      // Charged species
      const BoutReal charge1 = get<BoutReal>(species1["charge"]) * SI::qe; // in Coulombs
      
      // Copy the iterator, so we don't iterate over the
      // lower half of the matrix, but start at the diagonal
      for (std::map<std::string, Options>::const_iterator kv2 = kv1;
           kv2 != std::end(children); ++kv2) {
        if (kv2->first == "e")
          continue; // Skip electrons
        
        Options& species2 = allspecies[kv2->first];
        
        // Note: Here species1 could be equal to species2
        
        // If temperature isn't set, assume zero
        const Field3D temperature2 =
          species2.isSet("temperature") ? get<Field3D>(species2["temperature"]) : 0.0;
        const BoutReal mass2 = get<BoutReal>(species2["AA"]) * SI::Mp; // in Kg
        
        if (species2.isSet("charge")) {
          // Both charged species
          const BoutReal charge2 = get<BoutReal>(species2["charge"]) * SI::qe; // in Coulombs

          
        } else {
          // species1 charged, species2 neutral
          
        }
      }
    } else {
      // species1 neutral

      // Copy the iterator, so we don't iterate over the
      // lower half of the matrix, but start at the diagonal
      for (std::map<std::string, Options>::const_iterator kv2 = kv1;
           kv2 != std::end(children); ++kv2) {
        if (kv2->first == "e")
          continue; // Skip electrons
        
        Options& species2 = allspecies[kv2->first];
        
        // Note: Here species1 could be equal to species2
        
        // If temperature isn't set, assume zero
        const Field3D temperature2 =
          species2.isSet("temperature") ? get<Field3D>(species2["temperature"]) : 0.0;
        const BoutReal mass2 = get<BoutReal>(species2["AA"]) * SI::Mp; // in Kg
        
        if (species2.isSet("charge")) {
          // species1 neutral, species2 charged
          const BoutReal charge2 = get<BoutReal>(species2["charge"]) * SI::qe; // in Coulombs
          
          
        } else {
          // Both species neutral
          
        }
      }
      
    }     
  }
}

