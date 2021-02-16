#pragma once
#ifndef NEUTRAL_PARALLEL_DIFFUSION_H
#define NEUTRAL_PARALLEL_DIFFUSION_H

#include "component.hxx"

/// Hydrogen charge exchange total rate coefficient
///
/// Reaction 0.1T from Amjuel (p38)
///
/// Scaled to different isotope masses and finite neutral particle
/// temperatures by using the effective temperature (Amjuel p43)
///
/// T_eff = (M/M_1)T_1 + (M/M_2)T_2
///
struct HydrogenChargeExchange : public Component {
  ///
  /// @param alloptions Settings, which should include:
  ///        - units
  ///          - eV
  ///          - inv_meters_cubed
  ///          - seconds
  HydrogenChargeExchange(std::string, Options& alloptions, Solver*) {
    // Get the units
    const auto& units = alloptions["units"];
    Tnorm = get<BoutReal>(units["eV"]);
    Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
    FreqNorm = 1. / get<BoutReal>(units["seconds"]);
  }

protected:
  BoutReal Tnorm, Nnorm, FreqNorm; ///< Normalisations

  /// Calculate the charge exchange cross-section
  ///
  /// atom1 + ion1 -> atom2 + ion2
  ///
  /// and transfer of mass, momentum and energy from:
  ///
  /// atom1 -> ion2, ion1 -> atom2
  ///
  /// Assumes that both atom1 and ion1 have:
  ///   - AA
  ///   - density
  ///   - velocity
  ///   - temperature
  ///
  /// Sets in all species:
  ///   - density_source     [If atom1 != atom2 or ion1 != ion2]
  ///   - momentum_source
  ///   - energy_source
  ///
  void calculate_rates(Options& atom1, Options& ion1, Options& atom2, Options& ion2);
};

/// Hydrogen charge exchange
/// Templated on a char to allow 'h', 'd' and 't' species to be treated with the same code
///
/// @tparam Isotope1   The isotope ('h', 'd' or 't') of the initial atom
/// @tparam Isotope2   The isotope ('h', 'd' or 't') of the initial ion
template <char Isotope1, char Isotope2>
struct HydrogenChargeExchangeIsotope : public HydrogenChargeExchange {
  HydrogenChargeExchangeIsotope(std::string name, Options& alloptions, Solver* solver)
      : HydrogenChargeExchange(name, alloptions, solver) {}

  void transform(Options& state) override {
    calculate_rates(state["species"][{Isotope1}],       // e.g. "h"
                    state["species"][{Isotope2, '+'}],  // e.g. "d+"
                    state["species"][{Isotope2}],       // e.g. "d"
                    state["species"][{Isotope1, '+'}]); // e.g. "h+"
  }
};

namespace {
/// Register three components, one for each hydrogen isotope
/// so no isotope dependence included.
RegisterComponent<HydrogenChargeExchangeIsotope<'h', 'h'>>
    register_cx_hh("h + h+ -> h+ + h");
RegisterComponent<HydrogenChargeExchangeIsotope<'d', 'd'>>
    register_cx_dd("d + d+ -> d+ + d");
RegisterComponent<HydrogenChargeExchangeIsotope<'t', 't'>>
    register_cx_tt("t + t+ -> t+ + t");

// Charge exchange between different isotopes
RegisterComponent<HydrogenChargeExchangeIsotope<'h', 'd'>>
    register_cx_hd("h + d+ -> h+ + d");
RegisterComponent<HydrogenChargeExchangeIsotope<'d', 'h'>>
    register_cx_dh("d + h+ -> d+ + h");

RegisterComponent<HydrogenChargeExchangeIsotope<'h', 't'>>
    register_cx_ht("h + t+ -> h+ + t");
RegisterComponent<HydrogenChargeExchangeIsotope<'t', 'h'>>
    register_cx_th("t + h+ -> t+ + h");

RegisterComponent<HydrogenChargeExchangeIsotope<'d', 't'>>
    register_cx_dt("d + t+ -> d+ + t");
RegisterComponent<HydrogenChargeExchangeIsotope<'t', 'd'>>
    register_cx_td("t + d+ -> t+ + d");
} // namespace

#endif
