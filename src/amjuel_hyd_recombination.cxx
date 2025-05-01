#include "../include/amjuel_hyd_recombination.hxx"

void AmjuelHydRecombination::calculate_rates(
    Options& electron, Options& atom, Options& ion, Field3D& reaction_rate,
    Field3D& momentum_exchange, Field3D& energy_exchange, Field3D& energy_loss,
    BoutReal& rate_multiplier, BoutReal& radiation_multiplier) {
  electron_reaction(electron, ion, atom, get_rate_coeffs(), get_rad_coeffs(),
                    get_electron_heating(), 
                    reaction_rate, momentum_exchange, energy_exchange, energy_loss,
                    rate_multiplier, radiation_multiplier);
}
