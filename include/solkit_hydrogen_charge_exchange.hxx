#pragma once
#ifndef SOLKIT_HYDROGEN_CHARGE_EXCHANGE_H
#define SOLKIT_HYDROGEN_CHARGE_EXCHANGE_H

#include "component.hxx"
#include "hermes_utils.hxx"
#include "reaction.hxx"
/// SOL-KiT Hydrogen charge exchange total rate coefficient
///
struct SOLKITHydrogenChargeExchange : public hermes::ReactionBase {
  ///
  /// @param alloptions Settings, which should include:
  ///        - units
  ///          - inv_meters_cubed
  ///          - seconds
  SOLKITHydrogenChargeExchange(std::string name, Options& alloptions, Solver*)
      : hermes::ReactionBase(std::move(name),
                             {readIfSet("species:{sp}:velocity"),
                              readOnly("species:{sp}:{readvals}"),
                              readWrite("species:{sp}:momentum_source")}) {
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
  void calculate_rates(GuardedOptions atom, GuardedOptions ion);

protected:
  BoutReal Nnorm, rho_s0; ///< Normalisations
};

/// Hydrogen charge exchange
/// Templated on a char to allow 'h', 'd' and 't' species to be treated with the same code
///
/// @tparam Isotope   The isotope ('h', 'd' or 't') of the atom and ion
template <char Isotope>
struct SOLKITHydrogenChargeExchangeIsotope : public SOLKITHydrogenChargeExchange {
  SOLKITHydrogenChargeExchangeIsotope(std::string name, Options& alloptions,
                                      Solver* solver)
      : SOLKITHydrogenChargeExchange(name, alloptions, solver) {
    substitutePermissions("sp", {{Isotope}, {Isotope, '+'}});
    substitutePermissions("readvals", {"AA", "density"});
  }

  //FIXME: Need to figure out how to consturct type name at compile-time
  static constexpr char type[] = {'s', 'o',     'l', 'k',     'i', 't', ' ',     Isotope,
                                  ' ', '+',     ' ', Isotope, '+', ' ', '-',     '>',
                                  ' ', Isotope, '+', ' ',     '+', ' ', Isotope, '\0'};

  std::string typeName() const final { return type; }

private:
  void transform_impl(GuardedOptions& state) override {
    calculate_rates(state["species"][{Isotope}],       // e.g. "h"
                    state["species"][{Isotope, '+'}]); // e.g. "d+"
  }
};

namespace {
/// Register three components, one for each hydrogen isotope
RegisterComponent<SOLKITHydrogenChargeExchangeIsotope<'h'>> register_solkit_cx_hh;
RegisterComponent<SOLKITHydrogenChargeExchangeIsotope<'d'>> register_solkit_cx_dd;
RegisterComponent<SOLKITHydrogenChargeExchangeIsotope<'t'>> register_solkit_cx_tt;
} // namespace

#endif // SOLKIT_HYDROGEN_CHARGE_EXCHANGE_H
