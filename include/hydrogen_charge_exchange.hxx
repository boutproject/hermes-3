#pragma once
#ifndef HYDROGEN_CHARGE_EXCHANGE_H
#define HYDROGEN_CHARGE_EXCHANGE_H

#include <bout/constants.hxx>

#include "amjuel_reaction.hxx"
#include "component.hxx"

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
  HydrogenChargeExchange(std::string name, Options& alloptions, Solver*)
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
    std::string ion_reactant{Isotope2, '+'}; // rename to ion_out

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
        /*
         Simpler case of same isotopes  - No net particle source/sink; atoms lose
         * <atom_mom>, gain <ion_mom>
         */
      } else {
        // Different isotope => particle source, second momentum & energy channel

        // Need F2 = -ion_mom - CHECK
        add_diagnostic(
            atom_reactant, fmt::format("F{:s}{:s}+_cx", ion_product, atom_reactant),
            fmt::format("Momentum transfer to {:s} from {:s} due to CX with {:s}",
                        ion_product, atom_product, atom_reactant),
            ReactionDiagnosticType::momentum_src, "hydrogen_charge_exchange", identity,
            "momentum transfer");

        // Need E2 = -ion_energy - CHECK
        add_diagnostic(
            atom_reactant, fmt::format("E{:s}{:s}+_cx", ion_product, atom_reactant),
            fmt::format("Energy transfer to {:s} from {:s} due to CX with {:s}",
                        ion_product, atom_product, atom_reactant),
            ReactionDiagnosticType::energy_src, "hydrogen_charge_exchange", identity,
            "energy transfer");

        add_diagnostic(
            atom_reactant, fmt::format("S{:s}{:s}+_cx", atom_reactant, ion_product),
            fmt::format("Particle transfer to {:s} from {:s} due to CX with {:s}",
                        atom_reactant, ion_reactant, ion_product),
            ReactionDiagnosticType::density_src, "hydrogen_charge_exchange", identity,
            "particle transfer");
      }

      // Always add atom1 momentum source
      // Need transform such that
      // F = ion_mom - atom_mom for symmetric
      // F = -atom_mom for non-symmetric
      add_diagnostic(
          atom_reactant, fmt::format("F{:s}{:s}+_cx", atom_reactant, ion_product),
          fmt::format("Momentum transfer to {:s} from {:s} due to CX with {:s}",
                      atom_reactant, ion_reactant, ion_product),
          ReactionDiagnosticType::momentum_src, "hydrogen_charge_exchange",
          default_transformer);

      // Always add atom1 energy source
      // Need transform such that
      // E = ion_energy - atom_energy for symmetric
      // E = -atom_energy for non-symmetric
      add_diagnostic(atom_reactant,
                     fmt::format("E{:s}{:s}+_cx", atom_reactant, ion_product),
                     fmt::format("Energy transfer to {:s} from {:s} due to CX with {:s}",
                                 atom_reactant, ion_reactant, ion_product),
                     ReactionDiagnosticType::energy_src, "hydrogen_charge_exchange",
                     default_transformer, "energy transfer");

      // Always add CX collision frequency
      // Need transform such that
      // K = atom_rate
      add_diagnostic(atom_reactant,
                     fmt::format("K{:s}{:s}+_cx", atom_reactant, ion_product),
                     fmt::format("Collision frequency of CX of {:s} and {:s} producing "
                                 "{:s} and {:s}. Note Kab != Kba",
                                 atom_reactant, ion_reactant, ion_product, atom_product),
                     ReactionDiagnosticType::collision_freq, "hydrogen_charge_exchange",
                     default_transformer);
    }
  }

protected:
  virtual void transform_additional(Options& state, Field3D& reaction_rate) override {
    std::vector<std::string> ion_reactant =
        this->parser->get_species(species_filter::reactants, species_filter::ion);
    ASSERT1(ion_reactant.size() == 1);
    std::vector<std::string> ion_product =
        this->parser->get_species(species_filter::products, species_filter::ion);
    ASSERT1(ion_product.size() == 1);
    std::vector<std::string> heavy_reactant_species =
        this->parser->get_species(species_filter::reactants, species_filter::heavy);

    std::vector<std::string> neutral_reactant =
        this->parser->get_species(species_filter::reactants, species_filter::neutral);
    ASSERT1(neutral_reactant.size() == 1);
    std::vector<std::string> neutral_product =
        this->parser->get_species(species_filter::products, species_filter::neutral);
    ASSERT1(neutral_product.size() == 1);

    Options& atom1 = state["species"][neutral_reactant[0]];
    Options& ion1 = state["species"][ion_reactant[0]];
    Options& atom2 = state["species"][neutral_product[0]];
    Options& ion2 = state["species"][ion_product[0]];

    // Masses of initial atom and ion
    const BoutReal Aatom = get<BoutReal>(atom1["AA"]);
    // ASSERT1(get<BoutReal>(ion2["AA"]) == Aatom); // Check that the mass is consistent

    const BoutReal Aion = get<BoutReal>(ion1["AA"]);
    // ASSERT1(get<BoutReal>(atom2["AA"]) == Aion); // Check that the mass is consistent

    // Temporarily calculate sigmav here. Stop gap until collision_frequencies are
    // calculated in Reaction::transform
    Field3D Teff;
    calc_Teff(state, heavy_reactant_species, Teff);
    const Field3D lnT;
    Field3D sigmav(emptyFrom(Teff));
    BOUT_FOR(i, Teff.getRegion("RGN_NOBNDRY")) {
      sigmav[i] =
          this->eval_sigma_v_T(Teff[i]) * (1e-6 * Nnorm / FreqNorm) * rate_multiplier;
    }

    const Field3D Natom = floor(get<Field3D>(atom1["density"]), 1e-5);
    const Field3D Nion = floor(get<Field3D>(ion1["density"]), 1e-5);

    auto atom1_velocity = get<Field3D>(atom1["velocity"]);
    auto ion1_velocity = get<Field3D>(ion1["velocity"]);

    // Frictional heating: Friction force between ions and atoms
    // converts kinetic energy to thermal energy
    //
    // This handles the general case that ion1 != ion2
    // and atom1 != atom2

    auto ion2_velocity = get<Field3D>(ion2["velocity"]);
    add(ion2["energy_source"],
        0.5 * Aatom * reaction_rate * SQ(ion2_velocity - atom1_velocity));

    auto atom2_velocity = get<Field3D>(atom2["velocity"]);
    add(atom2["energy_source"],
        0.5 * Aion * reaction_rate * SQ(atom2_velocity - ion1_velocity));

    // Update collision frequency for the two colliding species
    Field3D atom_rate = Nion * sigmav; // [s^-1]
    Field3D ion_rate = Natom * sigmav; // [s^-1]

    // Add to total collision frequency
    add(atom1["collision_frequency"], atom_rate);
    add(ion1["collision_frequency"], ion_rate);

    // Set individual collision frequencies
    set(atom1["collision_frequencies"]
             [atom1.name() + std::string("_") + ion1.name() + std::string("_cx")],
        atom_rate);
    set(ion1["collision_frequencies"]
            [ion1.name() + std::string("_") + atom1.name() + std::string("_cx")],
        ion_rate);
  }

private:
  bool no_neutral_cx_mom_gain; ///< Make CX behave as in diffusive neutrals?
};

// void outputVars(Options& state) override {
//   AUTO_TRACE();
//   // Normalisations
//   auto Nnorm = get<BoutReal>(state["Nnorm"]);
//   auto Tnorm = get<BoutReal>(state["Tnorm"]);
//   BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation
//   auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
//   auto Cs0 = get<BoutReal>(state["Cs0"]);

//   if (diagnose) {
//     // Save particle, momentum and energy channels

//     std::string atom1{Isotope1};
//     std::string ion1{Isotope1, '+'};
//     std::string atom2{Isotope2};
//     std::string ion2{Isotope2, '+'};

//     set_with_attrs(state[{'F', Isotope1, Isotope2, '+', '_', 'c', 'x'}], // e.g Fhd+_cx
//                    F,
//                    {{"time_dimension", "t"},
//                     {"units", "kg m^-2 s^-2"},
//                     {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
//                     {"standard_name", "momentum transfer"},
//                     {"long_name", (std::string("Momentum transfer to ") + atom1
//                                    + " from " + ion1 + " due to CX with " + ion2)},
//                     {"source", "hydrogen_charge_exchange"}});

//     set_with_attrs(state[{'E', Isotope1, Isotope2, '+', '_', 'c', 'x'}], // e.g Edt+_cx
//                    E,
//                    {{"time_dimension", "t"},
//                     {"units", "W / m^3"},
//                     {"conversion", Pnorm * Omega_ci},
//                     {"standard_name", "energy transfer"},
//                     {"long_name", (std::string("Energy transfer to ") + atom1 + " from
//                     "
//                                    + ion1 + " due to CX with " + ion2)},
//                     {"source", "hydrogen_charge_exchange"}});

//     set_with_attrs(state[{'K', Isotope1, Isotope2, '+', '_', 'c', 'x'}], // e.g Kdt+_cx
//                    atom_rate,
//                    {{"time_dimension", "t"},
//                     {"units", "s^-1"},
//                     {"conversion", Omega_ci},
//                     {"standard_name", "collision frequency"},
//                     {"long_name", (std::string("Collision frequency of CX of ") + atom1
//                                    + " and " + ion1 + " producing " + ion2 + " and "
//                                    + atom2 + ". Note Kab != Kba")},
//                     {"source", "hydrogen_charge_exchange"}});

//     if (Isotope1 != Isotope2) {
//       // Different isotope => particle source, second momentum & energy channel
//       set_with_attrs(
//           state[{'F', Isotope2, '+', Isotope1, '_', 'c', 'x'}], // e.g Fd+h_cx
//           F2,
//           {{"time_dimension", "t"},
//            {"units", "kg m^-2 s^-2"},
//            {"conversion", SI::Mp * Nnorm * Cs0 * Omega_ci},
//            {"standard_name", "momentum transfer"},
//            {"long_name", (std::string("Momentum transfer to ") + ion2 + " from " +
//            atom2
//                           + " due to CX with " + atom1)},
//            {"source", "hydrogen_charge_exchange"}});

//       set_with_attrs(
//           state[{'E', Isotope2, '+', Isotope1, '_', 'c', 'x'}], // e.g Et+d_cx
//           E2,
//           {{"time_dimension", "t"},
//            {"units", "W / m^3"},
//            {"conversion", Pnorm * Omega_ci},
//            {"standard_name", "energy transfer"},
//            {"long_name", (std::string("Energy transfer to ") + ion2 + " from " + atom2
//                           + " due to CX with " + atom1)},
//            {"source", "hydrogen_charge_exchange"}});

//   }
// }
// };

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
