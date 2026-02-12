#pragma once
#ifndef HYDROGEN_CHARGE_EXCHANGE_H
#define HYDROGEN_CHARGE_EXCHANGE_H

#include <bout/constants.hxx>

#include "amjuel_reaction.hxx"
#include "component.hxx"
#include "reaction.hxx"

/**
* @brief Reaction component to handle Hydrogen isotope charge exchange.
*        p + H(1s) -> H(1s) + p
*        Templated on a char to allow 'h', 'd' and 't' to be treated with the same code
*
* @warning If this reaction is included then ion_neutral collisions should probably be
disabled in the `collisions` component, to avoid double-counting.
*
* @tparam Isotope1 The isotope ('h', 'd' or 't') of the reactant atom
* @tparam Isotope2 The isotope ('h', 'd' or 't') of the reactant ion
*
* i.e.
* atomR    +   ionR    ->   ionP     +    atomP
* Isotope1 + Isotope2+ -> Isotope1+  +  Isotope2
*
* @details
* Calculate the charge exchange cross-section for a reaction
*   atomR + ionR -> atomP + ionP
*   and handle (via Reaction::transform_impl) transfer of mass, momentum and energy from:
*   atomR -> ionP, ionR -> atomP
*
* Assumes that both atomR and ionR have:
*   - AA
*   - density
*   - velocity
*   - temperature
*
* Rate coefficients are computed from Reaction 3.1.8 from Amjuel (p43) using an effective
* temperature; see RateHelper::calc_Teff.
*
* Sets in all species (via Reaction::transform_impl and transform_additional):
*   - density_source     [If atomR != atomP or ionR != ionP]
*   - momentum_source
*   - energy_source
*
* Most of the source terms are handled by the Reaction base class.
* The transform_additional method in this class is used to add frictional heating to the
* energy source. transform_additional also sets collision_frequencies for atomR and ionR,
* which are calculated in Reaction::transform via the RateHelper class.
*
* Diagnostics
* -----------
*
* If diagnose = true is set in the options, then the following diagnostics are saved:
*   - F<Isotope1><Isotope2>+_cx  (e.g. Fhd+_cx) the momentum added to Isotope1 atoms due
*                                due to charge exchange with Isotope2 ions.
*                                There is a corresponding loss of momentum for the
*                                Isotope1 ions d/dt(NVh)  = ... + Fhd+_cx   // Atom
*                                momentum source d/dt(NVh+) = ... - Fhd+_cx   // Ion
*                                momentum sink
*   - E<Isotope1><Isotope2>+_cx  Energy added to Isotope1 atoms due to charge exchange
*   with
*                                Isotope2 ions. This contributes to two pressure
*                                equations d/dt(3/2 Ph)  = ... + Ehd+_cx d/dt(3/2 Ph+) =
*                                ... - Ehd+_cx
*
* If Isotope1 != Isotope2 then there is also the source of energy for Isotope2 atoms
* and a source of particles:
*   - F<Isotope2>+<Isotope1>_cx  Source of momentum for Isotope2 ions, sink for Isotope2
*   atoms
*   - E<Isotope2>+<Isotope1>_cx  Source of energy for Isotope2 ions, sink for Isotope2
*   atoms
*   - S<Isotope1><Isotope2>+_cx  Source of Isotope1 atoms due to charge exchange with
*   Isotope2 ions
*                                Note: S<Isotope2><Isotope1>+_cx =
*                                -S<Isotope1><Isotope2>+_cx For example Shd+_cx
*                                contributes to four density equations: d/dt(Nh)  = ...
*                                + Shd+_cx d/dt(Nh+) = ... - Shd+_cx d/dt(Nd)  = ... -
*                                Shd+_cx d/dt(Nd+) = ... + Shd+_cx
*/

template <char Isotope1, char Isotope2>
struct HydrogenChargeExchange : public AmjuelReaction {
  HydrogenChargeExchange([[maybe_unused]] std::string name, Options& alloptions, Solver*)
      : AmjuelReaction(name, "cx", "H.2_3.1.8", alloptions) {

    this->do_parallel_averaging = false;
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
    std::string atomR{Isotope1};
    std::string ionP{Isotope1, '+'};
    std::string atomP{Isotope2};
    std::string ionR{Isotope2, '+'};

    // Incoming neutral energy => outgoing ion
    this->set_energy_channel_weight(atomR, ionP, 1.0);
    this->set_energy_channel_weight(atomR, atomP, 0.0);
    // Incoming ion energy => outgoing neutral
    this->set_energy_channel_weight(ionR, ionP, 0.0);
    this->set_energy_channel_weight(ionR, atomP, 1.0);

    // Incoming neutral momentum => outgoing ion
    this->set_momentum_channel_weight(atomR, ionP, 1.0);
    this->set_momentum_channel_weight(atomR, atomP, 0.0);
    // Incoming ion momentum => outgoing neutral
    this->set_momentum_channel_weight(ionR, ionP, 0.0);
    if (this->no_neutral_cx_mom_gain) {
      this->set_momentum_channel_weight(ionR, atomP, 0.0);
    } else {
      this->set_momentum_channel_weight(ionR, atomP, 1.0);
    }

    // Set permissions
    setPermissions(readOnly("species:{reactant}:{react_vals}"));
    setPermissions(readOnly("species:{sp}:{read_vals}"));
    setPermissions(readWrite("species:{sp}:{writevals}"));
    setPermissions(readWrite("species:{atom}:collision_frequencies:{atom}_{ion}_cx"));
    setPermissions(readWrite("species:{ion}:collision_frequencies:{ion}_{atom}_cx"));

    std::vector<std::string> writevals = {"momentum_source", "energy_source"};
    if constexpr (Isotope1 != Isotope2) {
      writevals.push_back("density_source");
    }
    substitutePermissions("reactant", {{Isotope1}, {Isotope2, '+'}});
    substitutePermissions("react_vals", {"density", "temperature"});
    substitutePermissions("read_vals", {"AA", "velocity"});
    substitutePermissions("sp",
                          {{Isotope1}, {Isotope2, '+'}, {Isotope1, '+'}, {Isotope2}});
    substitutePermissions("writevals", writevals);
    substitutePermissions("atom", {{Isotope1}});
    substitutePermissions("ion", {{Isotope2, '+'}});

    if (diagnose) {
      // Set up diagnostics to be written to dump file
      DiagnosticTransformerType default_transformer = identity;

      if constexpr (Isotope1 != Isotope2) {
        // Different isotope => particle source, second momentum & energy channel
        add_diagnostic(
            atomR, fmt::format("S{:s}{:s}_cx", atomR, ionR),
            fmt::format("Particle transfer to {:s} from {:s} due to CX with {:s}", atomR,
                        ionP, ionR),
            ReactionDiagnosticType::density_src, "hydrogen_charge_exchange", identity,
            "particle transfer");

        // F2 = -ion-mom (key by atom product)
        add_diagnostic(
            atomP, fmt::format("F{:s}{:s}_cx", ionR, atomR),
            fmt::format("Momentum transfer to {:s} from {:s} due to CX with {:s}", ionR,
                        atomP, atomR),
            ReactionDiagnosticType::momentum_src, "hydrogen_charge_exchange", negate,
            "momentum transfer");
        // E2 = -ion-energy (key by atom product)
        add_diagnostic(
            atomP, fmt::format("E{:s}{:s}_cx", ionR, atomR),
            fmt::format("Energy transfer to {:s} from {:s} due to CX with {:s}", ionR,
                        atomP, atomR),
            ReactionDiagnosticType::energy_src, "hydrogen_charge_exchange", negate,
            "energy transfer");

        // F = -atom_mom for non-symmetric  (key by ion product)
        add_diagnostic(
            ionP, fmt::format("F{:s}{:s}_cx", atomR, ionR),
            fmt::format("Momentum transfer to {:s} from {:s} due to CX with {:s}", atomR,
                        ionP, ionR),
            ReactionDiagnosticType::momentum_src, "hydrogen_charge_exchange", negate);

        // E = -atom_energy for non-symmetric  (key by ion product)
        add_diagnostic(
            ionP, fmt::format("E{:s}{:s}_cx", atomR, ionR),
            fmt::format("Energy transfer to {:s} from {:s} due to CX with {:s}", atomR,
                        ionP, ionR),
            ReactionDiagnosticType::energy_src, "hydrogen_charge_exchange", negate,
            "energy transfer");

      } else {

        // F = ion_mom - atom_mom for symmetric
        add_diagnostic(
            ionP, fmt::format("F{:s}{:s}_cx", atomR, ionR),
            fmt::format("Momentum transfer to {:s} from {:s} due to CX with {:s}", atomR,
                        ionP, ionR),
            ReactionDiagnosticType::momentum_src, "hydrogen_charge_exchange", negate);

        // E = ion_energy - atom_energy for symmetric
        add_diagnostic(
            ionP, fmt::format("E{:s}{:s}_cx", atomR, ionR),
            fmt::format("Energy transfer to {:s} from {:s} due to CX with {:s}", atomR,
                        ionP, ionR),
            ReactionDiagnosticType::energy_src, "hydrogen_charge_exchange", negate,
            "energy transfer");
      }

      // Always add CX collision frequency
      add_diagnostic(atomR, fmt::format("K{:s}{:s}_cx", atomR, ionR),
          fmt::format("Collision frequency of CX of {:s} and {:s} producing "
                      "{:s} and {:s}. Note Kab != Kba",
                                 atomR, ionP, ionR, atomP),
                     ReactionDiagnosticType::collision_freq, "hydrogen_charge_exchange",
                     identity);
    }
  }

protected:
  void transform_additional(GuardedOptions& state, const RateData& rate_data) override {
    std::string ionR_name =
        this->parser->get_single_species(species_filter::reactants, species_filter::ion);
    std::string ionP_name =
        this->parser->get_single_species(species_filter::products, species_filter::ion);
    std::string atomR_name = this->parser->get_single_species(species_filter::reactants,
                                                              species_filter::neutral);
    std::string atomP_name = this->parser->get_single_species(species_filter::products,
                                                              species_filter::neutral);

    GuardedOptions atomR = state["species"][atomR_name];
    GuardedOptions ionR = state["species"][ionR_name];
    GuardedOptions atomP = state["species"][atomP_name];
    GuardedOptions ionP = state["species"][ionP_name];

    // Masses, veloci of initial atom and ion

    auto atomR_velocity = get<Field3D>(atomR["velocity"]);

    // Frictional heating: Friction force between ions and atoms
    // converts kinetic energy to thermal energy
    //
    // This handles the general case that ionR != ionP
    // and atomR != atomP
    const BoutReal AatomR = get<BoutReal>(atomR["AA"]);
    auto ionP_velocity = get<Field3D>(ionP["velocity"]);
    add(ionP["energy_source"],
        0.5 * AatomR * rate_data.rate * SQ(ionP_velocity - atomR_velocity));

    const BoutReal AionR = get<BoutReal>(ionR["AA"]);
    auto ionR_velocity = get<Field3D>(ionR["velocity"]);
    auto atom_p_velocity = get<Field3D>(atomP["velocity"]);
    add(atomP["energy_source"],
        0.5 * AionR * rate_data.rate * SQ(atom_p_velocity - ionR_velocity));

    // Set 'collision_frequencies' values
    std::string atom_coll_freq_key =
        fmt::format("collision_frequencies:{:s}_{:s}_cx", atomR_name, ionR_name);
    update_source<set>(state, atomR_name, ReactionDiagnosticType::collision_freq,
                       atom_coll_freq_key, rate_data.coll_freq(atomR_name));

    std::string ion_coll_freq_key =
        fmt::format("collision_frequencies:{:s}_{:s}_cx", ionR_name, atomR_name);
    update_source<set>(state, ionR_name, ReactionDiagnosticType::collision_freq,
                       ion_coll_freq_key, rate_data.coll_freq(ionR_name));
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
