#pragma once
#ifndef HYDROGEN_CHARGE_EXCHANGE_H
#define HYDROGEN_CHARGE_EXCHANGE_H

#include <bout/constants.hxx>

#include "amjuel_reaction.hxx"
#include "component.hxx"
#include "reaction.hxx"

/**
* @brief Reaction component to handle Hydrogen charge exchange.
*        p + H(1s) -> H(1s) + p
*        Templated on a char to allow 'h', 'd' and 't' to be treated with the same code
*
* @warning If this reaction is included then ion_neutral collisions should probably be
disabled in the `collisions` component, to avoid double-counting.

* @details Total rate coefficient is computed from: Reaction 3.1.8 from Amjuel (p43)
*
* Calculate the charge exchange cross-section for a reaction
*   atom1 + ion1 -> atom2 + ion2
*   and handle transfer of mass, momentum and energy from:
*   atom1 -> ion2, ion1 -> atom2
*
* Assumes that both atom1 and ion1 have:
*   - AA
*   - density
*   - velocity
*   - temperature
*
* Sets in all species:
*   - density_source     [If atom1 != atom2 or ion1 != ion2]
*   - momentum_source
*   - energy_source
*
* Modifies collision_frequency for atom1 and ion1
*
*  Diagnostic output (when diagnose = true is set in the options)
*   R            Reaction rate, transfer of particles in case of different isotopes
*   atom_mom     Momentum removed from atom1, added to ion2
*   ion_mom      Momentum removed from ion1, added to atom2
*   atom_energy  Energy removed from atom1, added to ion2
*   ion_energy   Energy removed from ion1, added to atom2
*
* @tparam Isotope1 The isotope ('h', 'd' or 't') of the reactant atom
* @tparam Isotope2 The isotope ('h', 'd' or 't') of the reactant ion
*
* i.e.
* atom     +   ion     ->   ion      +    atom
* Isotope1 + Isotope2+ -> Isotope1+  +  Isotope2
*/

// TODO: Add to docstring
///   - F<Isotope1><Isotope2>+_cx  (e.g. Fhd+_cx) the momentum added to Isotope1 atoms
///   due
///                                due to charge exchange with Isotope2 ions.
///                                There is a corresponding loss of momentum for the
///                                Isotope1 ions d/dt(NVh)  = ... + Fhd+_cx   // Atom
///                                momentum source d/dt(NVh+) = ... - Fhd+_cx   // Ion
///                                momentum sink
///   - E<Isotope1><Isotope2>+_cx  Energy added to Isotope1 atoms due to charge
///   exchange with
///                                Isotope2 ions. This contributes to two pressure
///                                equations d/dt(3/2 Ph)  = ... + Ehd+_cx d/dt(3/2
///                                Ph+)
///                                =
///                                ... - Ehd+_cx
///
/// If Isotope1 != Isotope2 then there is also the source of energy for Isotope2 atoms
/// and a source of particles:
///   - F<Isotope2>+<Isotope1>_cx  Source of momentum for Isotope2 ions, sink for
///   Isotope2 atoms
///   - E<Isotope2>+<Isotope1>_cx  Source of energy for Isotope2 ions, sink for
///   Isotope2 atoms
///   - S<Isotope1><Isotope2>+_cx  Source of Isotope1 atoms due to charge exchange
///   with Isotope2 ions
///                                Note: S<Isotope2><Isotope1>+_cx =
///                                -S<Isotope1><Isotope2>+_cx For example Shd+_cx
///                                contributes to four density equations: d/dt(Nh)  =
///                                ...
///                                + Shd+_cx d/dt(Nh+) = ... - Shd+_cx d/dt(Nd)  = ...
///                                - Shd+_cx d/dt(Nd+) = ... + Shd+_cx

template <char Isotope1, char Isotope2>
struct HydrogenChargeExchange : public AmjuelReaction {
  HydrogenChargeExchange([[maybe_unused]] std::string name, Options& alloptions, Solver*)
      : AmjuelReaction(name, "cx", "H.2_3.1.8", alloptions) {
    this->includes_sigma_v_e = false;
    /* This is useful for testing the impact of enabling the neutral momentum equation.
     * When set to true, CX behaves as if using diffusive neutrals but the neutral
     * transport still enjoys the full momentum equation treatment.
     */
    this->no_neutral_cx_mom_gain =
        alloptions[name]["no_neutral_cx_mom_gain"]
            .doc("If true, ion momentum in Hydrogen isotope CX is still lost, but not "
                 "given to the neutrals")
            .withDefault<bool>(false);

    // Options under neutral species of isotope 1 (on LHS of reaction)
    rate_multiplier = alloptions[{Isotope1}]["K_cx_multiplier"]
                          .doc("Scale the charge exchange rate by this factor")
                          .withDefault<BoutReal>(1.0);

    // Construct reactant, product names from template params
    std::string atom_reactant{Isotope1};
    std::string ion_product{Isotope1, '+'};
    std::string atom_product{Isotope2};
    std::string ion_reactant{Isotope2, '+'};

    // Incoming neutral energy => outgoing ion
    this->set_energy_channel_weight(atom_reactant, ion_product, 1.0);
    this->set_energy_channel_weight(atom_reactant, atom_product, 0.0);
    // Incoming ion energy => outgoing neutral
    this->set_energy_channel_weight(ion_reactant, ion_product, 0.0);
    this->set_energy_channel_weight(ion_reactant, atom_product, 1.0);

    // Incoming neutral momentum => outgoing ion
    this->set_momentum_channel_weight(atom_reactant, ion_product, 1.0);
    this->set_momentum_channel_weight(atom_reactant, atom_product, 0.0);
    // Incoming ion momentum => outgoing neutral
    this->set_momentum_channel_weight(ion_reactant, ion_product, 0.0);
    if (this->no_neutral_cx_mom_gain) {
      this->set_momentum_channel_weight(ion_reactant, atom_product, 0.0);
    } else {
      this->set_momentum_channel_weight(ion_reactant, atom_product, 1.0);
    }

    if (diagnose) {
      // Set up diagnostics to be written to dump file
      DiagnosticTransformerType default_transformer = identity;

      if constexpr (Isotope1 != Isotope2) {
        // Different isotope => particle source, second momentum & energy channel
        add_diagnostic(
            atom_reactant, fmt::format("S{:s}{:s}_cx", atom_reactant, ion_reactant),
            fmt::format("Particle transfer to {:s} from {:s} due to CX with {:s}",
                        atom_reactant, ion_product, ion_reactant),
            ReactionDiagnosticType::density_src, "hydrogen_charge_exchange", identity,
            "particle transfer");
      } else {
        /*
         Simpler case of same isotopes  - No net particle source/sink; atoms lose
         * <atom_mom>, gain <ion_mom>
         */
        // Need F2 = -ion_mom - CHECK
        add_diagnostic(
            atom_reactant, fmt::format("F{:s}{:s}_cx", ion_product, atom_reactant),
            fmt::format("Momentum transfer to {:s} from {:s} due to CX with {:s}",
                        ion_product, atom_product, atom_reactant),
            ReactionDiagnosticType::momentum_src, "hydrogen_charge_exchange", identity,
            "momentum transfer");

        // Need E2 = -ion_energy - CHECK
        add_diagnostic(
            atom_reactant, fmt::format("E{:s}{:s}_cx", ion_product, atom_reactant),
            fmt::format("Energy transfer to {:s} from {:s} due to CX with {:s}",
                        ion_product, atom_product, atom_reactant),
            ReactionDiagnosticType::energy_src, "hydrogen_charge_exchange", identity,
            "energy transfer");
      }

      // Always add atom1 momentum source
      // Need transform such that
      // F = ion_mom - atom_mom for symmetric
      // F = -atom_mom for non-symmetric
      add_diagnostic(
          atom_reactant, fmt::format("F{:s}{:s}_cx", atom_reactant, ion_product),
          fmt::format("Momentum transfer to {:s} from {:s} due to CX with {:s}",
                      atom_reactant, ion_reactant, ion_product),
          ReactionDiagnosticType::momentum_src, "hydrogen_charge_exchange",
          default_transformer);

      // Always add atom1 energy source
      // Need transform such that
      // E = ion_energy - atom_energy for symmetric
      // E = -atom_energy for non-symmetric
      add_diagnostic(atom_reactant,
                     fmt::format("E{:s}{:s}_cx", atom_reactant, ion_product),
                     fmt::format("Energy transfer to {:s} from {:s} due to CX with {:s}",
                                 atom_reactant, ion_reactant, ion_product),
                     ReactionDiagnosticType::energy_src, "hydrogen_charge_exchange",
                     default_transformer, "energy transfer");

      // Always add CX collision frequency
      // Need transform such that
      // K = atom_rate
      add_diagnostic(atom_reactant,
                     fmt::format("K{:s}{:s}_cx", atom_reactant, ion_product),
                     fmt::format("Collision frequency of CX of {:s} and {:s} producing "
                                 "{:s} and {:s}. Note Kab != Kba",
                                 atom_reactant, ion_reactant, ion_product, atom_product),
                     ReactionDiagnosticType::collision_freq, "hydrogen_charge_exchange",
                     default_transformer);
    }
  }

protected:
  void transform_additional(Options& state, RatesMap& rate_calc_results) override {
    std::string ion_reactant =
        this->parser->get_single_species(species_filter::reactants, species_filter::ion);
    std::string ion_product =
        this->parser->get_single_species(species_filter::products, species_filter::ion);
    std::string neutral_reactant = this->parser->get_single_species(
        species_filter::reactants, species_filter::neutral);
    std::string neutral_product = this->parser->get_single_species(
        species_filter::products, species_filter::neutral);

    Options& atom1 = state["species"][neutral_reactant];
    Options& ion1 = state["species"][ion_reactant];
    Options& atom2 = state["species"][neutral_product];
    Options& ion2 = state["species"][ion_product];

    // Masses of initial atom and ion
    const BoutReal Aatom = get<BoutReal>(atom1["AA"]);
    // ASSERT1(get<BoutReal>(ion2["AA"]) == Aatom); // Check that the mass is consistent

    const BoutReal Aion = get<BoutReal>(ion1["AA"]);
    // ASSERT1(get<BoutReal>(atom2["AA"]) == Aion); // Check that the mass is consistent

    auto atom1_velocity = get<Field3D>(atom1["velocity"]);
    auto ion1_velocity = get<Field3D>(ion1["velocity"]);

    // Frictional heating: Friction force between ions and atoms
    // converts kinetic energy to thermal energy
    //
    // This handles the general case that ion1 != ion2
    // and atom1 != atom2
    auto ion2_velocity = get<Field3D>(ion2["velocity"]);
    add(ion2["energy_source"],
        0.5 * Aatom * rate_calc_results["rate"] * SQ(ion2_velocity - atom1_velocity));

    auto atom2_velocity = get<Field3D>(atom2["velocity"]);
    add(atom2["energy_source"],
        0.5 * Aion * rate_calc_results["rate"] * SQ(atom2_velocity - ion1_velocity));

    // Set individual collision frequencies
    std::string neutral_coll_freq_key =
        fmt::format("{:s}:collision_frequency", neutral_reactant);
    set(atom1["collision_frequencies"]
             [atom1.name() + std::string("_") + ion1.name() + std::string("_cx")],
        rate_calc_results[neutral_coll_freq_key]);

    std::string ion_coll_freq_key = fmt::format("{:s}:collision_frequency", ion_reactant);
    set(ion1["collision_frequencies"]
            [ion1.name() + std::string("_") + atom1.name() + std::string("_cx")],
        rate_calc_results[ion_coll_freq_key]);
  }

private:
  bool no_neutral_cx_mom_gain; ///< Make CX behave as in diffusive neutrals?
};

namespace {
/// Register three components, one for each hydrogen isotope
/// so no isotope dependence included.
RegisterComponent<HydrogenChargeExchange<'h', 'h'>> register_cx_hh("h + h+ -> h+ + h");
RegisterComponent<HydrogenChargeExchange<'d', 'd'>> register_cx_dd("d + d+ -> d+ + d");
RegisterComponent<HydrogenChargeExchange<'t', 't'>> register_cx_tt("t + t+ -> t+ + t");

// Charge exchange between different isotopes
RegisterComponent<HydrogenChargeExchange<'h', 'd'>> register_cx_hd("h + d+ -> h+ + d");
RegisterComponent<HydrogenChargeExchange<'d', 'h'>> register_cx_dh("d + h+ -> d+ + h");

RegisterComponent<HydrogenChargeExchange<'h', 't'>> register_cx_ht("h + t+ -> h+ + t");
RegisterComponent<HydrogenChargeExchange<'t', 'h'>> register_cx_th("t + h+ -> t+ + h");

RegisterComponent<HydrogenChargeExchange<'d', 't'>> register_cx_dt("d + t+ -> d+ + t");
RegisterComponent<HydrogenChargeExchange<'t', 'd'>> register_cx_td("t + d+ -> t+ + d");
} // namespace

#endif // HYDROGEN_CHARGE_EXCHANGE_H
