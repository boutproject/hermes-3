#pragma once
#ifndef SOLKIT_HYDROGEN_CHARGE_EXCHANGE_H
#define SOLKIT_HYDROGEN_CHARGE_EXCHANGE_H

#include "component.hxx"

/// SOL-KiT Hydrogen charge exchange total rate coefficient
///
struct SOLKITHydrogenChargeExchange : public Component {
  ///
  /// @param alloptions Settings, which should include:
  ///        - units
  ///          - inv_meters_cubed
  ///          - seconds
  SOLKITHydrogenChargeExchange(std::string, Options& alloptions, Solver*) {
    // Get the units
    const auto& units = alloptions["units"];
    Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
    rho_s0 = get<BoutReal>(units["meters"]);
  }

  /// Calculate the charge exchange cross-section
  ///
  /// atom + ion -> atom + ion
  ///
  /// Assumes that both atom and ion have:
  ///   - AA
  ///   - density
  ///   - velocity
  ///
  /// Sets in all species:
  ///   - momentum_source
  ///
  void calculate_rates(Options& atom, Options& ion);

protected:
  BoutReal Nnorm, rho_s0; ///< Normalisations
};

/// Hydrogen charge exchange
/// Templated on a char to allow 'h', 'd' and 't' species to be treated with the same code
///
/// @tparam Isotope   The isotope ('h', 'd' or 't') of the atom and ion
template <char Isotope>
struct SOLKITHydrogenChargeExchangeIsotope : public SOLKITHydrogenChargeExchange {
  SOLKITHydrogenChargeExchangeIsotope(std::string name, Options& alloptions, Solver* solver)
      : SOLKITHydrogenChargeExchange(name, alloptions, solver) {}

  void transform(Options& state) override {
    calculate_rates(state["species"][{Isotope}],        // e.g. "h"
                    state["species"][{Isotope, '+'}]);  // e.g. "d+"
  }
};

namespace {
/// Register three components, one for each hydrogen isotope
RegisterComponent<SOLKITHydrogenChargeExchangeIsotope<'h'>>
    register_solkit_cx_hh("solkit h + h+ -> h+ + h");
RegisterComponent<SOLKITHydrogenChargeExchangeIsotope<'d'>>
    register_solkit_cx_dd("solkit d + d+ -> d+ + d");
RegisterComponent<SOLKITHydrogenChargeExchangeIsotope<'t'>>
    register_solkit_cx_tt("solkit t + t+ -> t+ + t");
} // namespace

#endif // SOLKIT_HYDROGEN_CHARGE_EXCHANGE_H
