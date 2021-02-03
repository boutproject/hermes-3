
#include "../include/ionisation.hxx"
#include "../include/integrate.hxx"

Ionisation::Ionisation(std::string name, Options &alloptions, Solver *) {

  // Get options for this component
  auto& options = alloptions[name];
  
  Eionize = options["Eionize"].doc("Ionisation energy cost").withDefault(30.);

  // Get the units
  const auto& units = alloptions["units"];
  Tnorm = get<BoutReal>(units["eV"]);
  Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
  FreqNorm = 1. / get<BoutReal>(units["seconds"]);

  // Normalise
  Eionize /= Tnorm;
}

void Ionisation::transform(Options &state) {
  // Get neutral atom properties
  Options& hydrogen = state["species"]["h"];
  Field3D Nn = get<Field3D>(hydrogen["density"]);
  Field3D Tn = get<Field3D>(hydrogen["temperature"]);
  Field3D Vn = get<Field3D>(hydrogen["velocity"]);
  auto AA = get<BoutReal>(hydrogen["AA"]);
  
  Options& electron = state["species"]["e"];
  Field3D Ne = get<Field3D>(electron["density"]);
  Field3D Te = get<Field3D>(electron["temperature"]);

  Options& ion = state["species"]["h+"];
  ASSERT1(AA == get<BoutReal>(ion["AA"]));

  Field3D reaction_rate = cellAverage(
      [&](BoutReal ne, BoutReal nn, BoutReal te) {
        return ne * nn * atomic_rates.ionisation(te * Tnorm) * Nnorm / FreqNorm;
      },
      Ne.getRegion("RGN_NOBNDRY"))(Ne, Nn, Te);

  // Particles move from hydrogen to ion
  subtract(hydrogen["density_source"],
           reaction_rate);

  add(ion["density_source"],
      reaction_rate);

  // Move momentum from hydrogen to ion
  Field3D momentum_exchange = reaction_rate * AA * Vn;
  
  subtract(hydrogen["momentum_source"],
           momentum_exchange);
  add(ion["momentum_source"],
      momentum_exchange);

  // Move energy from hydrogen to ion
  Field3D energy_exchange = reaction_rate * (3. / 2) * Tn;
  subtract(hydrogen["energy_source"],
           energy_exchange);
  add(ion["energy_source"],
      energy_exchange);

  // Radiate energy from electrons
  subtract(electron["energy_source"],
           Eionize * reaction_rate);
}
