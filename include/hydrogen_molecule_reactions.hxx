#pragma once
#ifndef HYDROGEN_MOLECULE_REACTIONS_H
#define HYDROGEN_MOLECULE_REACTIONS_H

#include "cx_reaction.hxx"
#include "molecular_reactions.hxx"
#include <string>

namespace hermes {

enum class HIsotope { h, d, t };

/**
 * @brief Reaction class for molecular charge exchange / ion conversion.
 *
 * @details Based on CXReaction.
 */
struct MolHCX : public CXReaction {
  /**
   * @brief Main constructor for MolHCX.
   *
   * @param name
   * @param options The options object
   */
  MolHCX(std::string name, Options& options);

  /**
   * @brief Constructor used by component factory.
   *
   * @param name
   * @param options The options object
   * @param solver The solver object for the simulation
   */
  MolHCX(std::string name, Options& options, Solver*);
};

/**
 * @brief Class to handle dissociation of molecular Hydrogen (and its isotopes).
 */
template <HIsotope isotope>
struct MolHDissociation : public Dissociation {
  /**
   * @brief Constructor for MolHDissociation.
   *
   * @param name
   * @param options The options object
   */
  MolHDissociation(std::string name, Options& options) : Dissociation(name, options) {}

  /**
   * @brief Constructor used by component factory.
   *
   * @param name
   * @param options The options object
   * @param solver The solver object for the simulation
   */
  MolHDissociation(std::string name, Options& options, Solver* solver)
      : Dissociation(name, options, solver) {}
};

/**
 * @brief Class to handle dissociative excitation of molecular Hydrogen (and its
 * isotopes).
 */
template <HIsotope isotope>
struct MolHDissociativeExc : public DissociativeExc {
  /**
   * @brief Constructor for MolHDissociativeExc.
   *
   * @param name
   * @param options The options object
   */
  MolHDissociativeExc(std::string name, Options& options)
      : DissociativeExc(name, options) {}

  /**
   * @brief Constructor used by component factory.
   *
   * @param name
   * @param options The options object
   * @param solver The solver object for the simulation
   */
  MolHDissociativeExc(std::string name, Options& options, Solver* solver)
      : DissociativeExc(name, options, solver) {}
};

/**
 * @brief Class to handle dissociative ionisation of molecular Hydrogen (and its
 * isotopes).
 */
template <HIsotope isotope>
struct MolHDissociativeIzn : public DissociativeIzn {
  /**
   * @brief Constructor for MolHDissociativeIzn.
   *
   * @param name
   * @param options The options object
   */
  MolHDissociativeIzn(std::string name, Options& options)
      : DissociativeIzn(name, options) {}

  /**
   * @brief Constructor used by component factory.
   *
   * @param name
   * @param options The options object
   * @param solver The solver object for the simulation
   */
  MolHDissociativeIzn(std::string name, Options& options, Solver* solver)
      : DissociativeIzn(name, options, solver) {}
};

/**
 * @brief Class to handle non-dissociative ionisation of molecular Hydrogen (and its
 * isotopes).
 */
template <HIsotope isotope>
struct MolHNonDissociativeIzn : public NonDissociativeIzn {
  /**
   * @brief Constructor for MolHNonDissociativeIzn.
   *
   * @param name
   * @param options The options object
   */
  MolHNonDissociativeIzn(std::string name, Options& options)
      : NonDissociativeIzn(name, options) {}

  /**
   * @brief Constructor used by component factory.
   *
   * @param name
   * @param options The options object
   * @param solver The solver object for the simulation
   */
  MolHNonDissociativeIzn(std::string name, Options& options, Solver* solver)
      : NonDissociativeIzn(name, options, solver) {}
};

/**
 * @brief Class to handle dissociative recombination of molecular Hydrogen (and its
 * isotopes).
 */
template <HIsotope isotope>
struct MolHDissociativeRec : public DissociativeRec {
  /**
   * @brief Constructor for MolHDissociativeRec.
   *
   * @param name
   * @param options The options object
   */
  MolHDissociativeRec(std::string name, Options& options)
      : DissociativeRec(name, options) {}

  /**
   * @brief Constructor used by component factory.
   *
   * @param name
   * @param options The options object
   * @param solver The solver object for the simulation
   */
  MolHDissociativeRec(std::string name, Options& options, Solver* solver)
      : DissociativeRec(name, options, solver) {}
};

} // namespace hermes

// Register molecular Hydrogen reactions
namespace {

// Non-dissociative ionisation of H isotope molecules
RegisterComponent<hermes::MolHNonDissociativeIzn<hermes::HIsotope::h>>
    register_h2_nondiss_izn_h("h2 + e -> h2+ + 2e");
RegisterComponent<hermes::MolHNonDissociativeIzn<hermes::HIsotope::d>>
    register_d2_nondiss_izn_d("d2 + e -> d2+ + 2e");
RegisterComponent<hermes::MolHNonDissociativeIzn<hermes::HIsotope::t>>
    register_t2_nondiss_izn_t("t2 + e -> t2+ + 2e");

// Dissociation of H isotope molecules
RegisterComponent<hermes::MolHDissociation<hermes::HIsotope::h>>
    register_h2_diss("h2 + e -> 2h + e");
RegisterComponent<hermes::MolHDissociation<hermes::HIsotope::d>>
    register_d2_diss("d2 + e -> 2d + e");
RegisterComponent<hermes::MolHDissociation<hermes::HIsotope::t>>
    register_t2_diss("t2 + e -> 2t + e");

// Dissociative ionisation of H isotope molecules
RegisterComponent<hermes::MolHDissociativeIzn<hermes::HIsotope::h>>
    register_h2_diss_izn_h("h2 + e -> h + h+ + 2e");
RegisterComponent<hermes::MolHDissociativeIzn<hermes::HIsotope::d>>
    register_d2_diss_izn_d("d2 + e -> d + d+ + 2e");
RegisterComponent<hermes::MolHDissociativeIzn<hermes::HIsotope::t>>
    register_t2_diss_izn_t("t2 + e -> t + t+ + 2e");

// Charge exchange between H isotope molecules and H isotope ions
RegisterComponent<hermes::MolHCX> register_h2_h_cx("h2 + h+ -> h2+ + h");
RegisterComponent<hermes::MolHCX> register_d2_d_cx("d2 + d+ -> d2+ + d");
RegisterComponent<hermes::MolHCX> register_t2_t_cx("t2 + t+ -> t2+ + t");

// Dissociative recombination of H isotope molecular ions
RegisterComponent<hermes::MolHDissociativeRec<hermes::HIsotope::h>>
    register_h2_plus_rec_h("h2+ + e -> 2h");
RegisterComponent<hermes::MolHDissociativeRec<hermes::HIsotope::d>>
    register_d2_plus_rec_d("d2+ + e -> 2d");
RegisterComponent<hermes::MolHDissociativeRec<hermes::HIsotope::t>>
    register_t2_plus_rec_t("t2+ + e -> 2t");

// Dissociative excitation of H isotope molecular ions
RegisterComponent<hermes::MolHDissociativeExc<hermes::HIsotope::h>>
    register_h2_plus_diss_exc_h("h2+ + e -> h + h+ + e");
RegisterComponent<hermes::MolHDissociativeExc<hermes::HIsotope::d>>
    register_d2_plus_diss_exc_d("d2+ + e -> d + d+ + e");
RegisterComponent<hermes::MolHDissociativeExc<hermes::HIsotope::t>>
    register_t2_plus_diss_exc_t("t2+ + e -> t + t+ + e");

// Dissociative ionisation of H isotope molecular ions
RegisterComponent<hermes::MolHDissociativeIzn<hermes::HIsotope::h>>
    register_h2_plus_diss_izn_h("h2+ + e -> 2h+ + 2e");
RegisterComponent<hermes::MolHDissociativeIzn<hermes::HIsotope::d>>
    register_d2_plus_diss_izn_d("d2+ + e -> 2d+ + 2e");
RegisterComponent<hermes::MolHDissociativeIzn<hermes::HIsotope::t>>
    register_t2_plus_diss_izn_t("t2+ + e -> 2t+ + 2e");

} // namespace

#endif // HYDROGEN_MOLECULE_REACTIONS_H
