#include "../include/amjuel_helium.hxx"

////////////////////////////////////////////////////////////
// e + he -> he+ + 2e

void AmjuelHeIonisation01::calculate_rates(Options& state, Field3D& reaction_rate,
                                           Field3D& momentum_exchange,
                                           Field3D& energy_exchange, Field3D& energy_loss,
                                           BoutReal& rate_multiplier,
                                           BoutReal& radiation_multiplier) {
  electron_reaction(state["species"]["e"],
                    state["species"]["he"],  // From helium atoms
                    state["species"]["he+"], // To helium ions
                    get_rate_coeffs(), get_rad_coeffs(),
                    get_electron_heating(), 
                    reaction_rate, momentum_exchange, energy_exchange, energy_loss,
                    rate_multiplier, radiation_multiplier

  );
}

////////////////////////////////////////////////////////////
// e + he+ -> he

void AmjuelHeRecombination10::calculate_rates(Options& state, Field3D& reaction_rate,
                                              Field3D& momentum_exchange,
                                              Field3D& energy_exchange,
                                              Field3D& energy_loss,
                                              BoutReal& rate_multiplier,
                                              BoutReal& radiation_multiplier) {
  electron_reaction(state["species"]["e"],
                    state["species"]["he+"], // From helium ions
                    state["species"]["he"],  // To helium atoms
                    get_rate_coeffs(), get_rad_coeffs(),
                    get_electron_heating(),
                    reaction_rate, momentum_exchange, energy_exchange, energy_loss,
                    rate_multiplier, radiation_multiplier);
}
