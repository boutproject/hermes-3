#pragma once
#ifndef HYDROGEN_MOLECULE_REACTIONS_H
#define HYDROGEN_MOLECULE_REACTIONS_H

#include "cx_reaction.hxx"
#include "elastic_collisions.hxx"
#include "molecular_reactions.hxx"
#include <string>

namespace hermes {

enum class HIsotope { h, d, t };

/**
 * @brief Class to handle charge exchange between Hydrogen (isotope) ions and molecules.
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

/**
 * @brief Templated hydrogen molecule elastic collision reaction.
 *
 * @details Child class of ElasticCollision, templated on hydrogen isotope.
 */
template <HIsotope isotope>
struct MolHElasticCollision : public ElasticCollision {
  /**
   * @brief Constructor for MolHElasticCollision.
   *
   * @param name The name of the reaction
   * @param options The options object
   */
  MolHElasticCollision(std::string name, Options& options)
      : ElasticCollision(name, options) {}

  /**
   * @brief Constructor used by component factory.
   *
   * @param name The name of the reaction
   * @param options The options object
   * @param solver The solver object for the simulation
   */
  MolHElasticCollision(std::string name, Options& options,
                       [[maybe_unused]] Solver* solver)
      : MolHElasticCollision(name, options) {}
};

} // namespace hermes

// Register molecular Hydrogen reactions
namespace {

// Non-dissociative ionisation of H isotope molecules
struct H2NonDissIzn : public hermes::MolHNonDissociativeIzn<hermes::HIsotope::h> {
  using MolHNonDissociativeIzn<hermes::HIsotope::h>::MolHNonDissociativeIzn;
  static constexpr auto type = "h2 + e -> h2+ + 2e";
};
RegisterComponent<H2NonDissIzn> register_h2_non_diss_izn;
struct D2NonDissIzn : public hermes::MolHNonDissociativeIzn<hermes::HIsotope::d> {
  using MolHNonDissociativeIzn<hermes::HIsotope::d>::MolHNonDissociativeIzn;
  static constexpr auto type = "d2 + e -> d2+ + 2e";
};
RegisterComponent<D2NonDissIzn> register_d2_non_diss_izn;
struct T2NonDissIzn : public hermes::MolHNonDissociativeIzn<hermes::HIsotope::t> {
  using MolHNonDissociativeIzn<hermes::HIsotope::t>::MolHNonDissociativeIzn;
  static constexpr auto type = "t2 + e -> t2+ + 2e";
};
RegisterComponent<T2NonDissIzn> register_t2_non_diss_izn;

// Dissociation of H isotope molecules
struct H2Dissociation : public hermes::MolHDissociation<hermes::HIsotope::h> {
  using MolHDissociation<hermes::HIsotope::h>::MolHDissociation;
  static constexpr auto type = "h2 + e -> 2h + e";
};
RegisterComponent<H2Dissociation> register_h2_dissociation;
struct D2Dissociation : public hermes::MolHDissociation<hermes::HIsotope::d> {
  using MolHDissociation<hermes::HIsotope::d>::MolHDissociation;
  static constexpr auto type = "d2 + e -> 2d + e";
};
RegisterComponent<D2Dissociation> register_d2_dissociation;
struct T2Dissociation : public hermes::MolHDissociation<hermes::HIsotope::t> {
  using MolHDissociation<hermes::HIsotope::t>::MolHDissociation;
  static constexpr auto type = "t2 + e -> 2t + e";
};
RegisterComponent<T2Dissociation> register_t2_dissociation;

// Dissociative ionisation of H isotope molecules
struct H2DissIzn : public hermes::MolHDissociativeIzn<hermes::HIsotope::h> {
  using MolHDissociativeIzn<hermes::HIsotope::h>::MolHDissociativeIzn;
  static constexpr auto type = "h2 + e -> h + h+ + 2e";
};
RegisterComponent<H2DissIzn> register_h2_diss_izn;
struct D2DissIzn : public hermes::MolHDissociativeIzn<hermes::HIsotope::d> {
  using MolHDissociativeIzn<hermes::HIsotope::d>::MolHDissociativeIzn;
  static constexpr auto type = "d2 + e -> d + d+ + 2e";
};
RegisterComponent<D2DissIzn> register_d2_diss_izn;
struct T2DissIzn : public hermes::MolHDissociativeIzn<hermes::HIsotope::t> {
  using MolHDissociativeIzn<hermes::HIsotope::t>::MolHDissociativeIzn;
  static constexpr auto type = "t2 + e -> t + t+ + 2e";
};
RegisterComponent<T2DissIzn> register_t2_diss_izn;

// Charge exchange between H isotope molecules and H isotope ions
struct H2CX : public hermes::MolHCX {
  using MolHCX::MolHCX;
  static constexpr auto type = "h2 + h+ -> h2+ + h";
};
RegisterComponent<H2CX> register_h2_cx;
struct D2CX : public hermes::MolHCX {
  using MolHCX::MolHCX;
  static constexpr auto type = "d2 + d+ -> d2+ + d";
};
RegisterComponent<D2CX> register_d2_cx;
struct T2CX : public hermes::MolHCX {
  using MolHCX::MolHCX;
  static constexpr auto type = "t2 + t+ -> t2+ + t";
};
RegisterComponent<T2CX> register_t2_cx;

// Dissociative recombination of H isotope molecular ions
struct H2PlusRec : public hermes::MolHDissociativeRec<hermes::HIsotope::h> {
  using MolHDissociativeRec<hermes::HIsotope::h>::MolHDissociativeRec;
  static constexpr auto type = "h2+ + e -> 2h";
};
RegisterComponent<H2PlusRec> register_h2_plus_rec;
struct D2PlusRec : public hermes::MolHDissociativeRec<hermes::HIsotope::d> {
  using MolHDissociativeRec<hermes::HIsotope::d>::MolHDissociativeRec;
  static constexpr auto type = "d2+ + e -> 2d";
};
RegisterComponent<D2PlusRec> register_d2_plus_rec;

struct T2PlusRec : public hermes::MolHDissociativeRec<hermes::HIsotope::t> {
  using MolHDissociativeRec<hermes::HIsotope::t>::MolHDissociativeRec;
  static constexpr auto type = "t2+ + e -> 2t";
};
RegisterComponent<T2PlusRec> register_t2_plus_rec;

// Dissociative excitation of H isotope molecular ions
struct H2PlusDissExc : public hermes::MolHDissociativeExc<hermes::HIsotope::h> {
  using MolHDissociativeExc<hermes::HIsotope::h>::MolHDissociativeExc;
  static constexpr auto type = "h2+ + e -> h + h+ + e";
};
RegisterComponent<H2PlusDissExc> register_h2_plus_diss_exc;
struct D2PlusDissExc : public hermes::MolHDissociativeExc<hermes::HIsotope::d> {
  using MolHDissociativeExc<hermes::HIsotope::d>::MolHDissociativeExc;
  static constexpr auto type = "d2+ + e -> d + d+ + e";
};
RegisterComponent<D2PlusDissExc> register_d2_plus_diss_exc;
struct T2PlusDissExc : public hermes::MolHDissociativeExc<hermes::HIsotope::t> {
  using MolHDissociativeExc<hermes::HIsotope::t>::MolHDissociativeExc;
  static constexpr auto type = "t2+ + e -> t + t+ + e";
};
RegisterComponent<T2PlusDissExc> register_t2_plus_diss_exc;

// Dissociative ionisation of H isotope molecular ions
struct H2PlusDissIzn : public hermes::MolHDissociativeIzn<hermes::HIsotope::h> {
  using MolHDissociativeIzn<hermes::HIsotope::h>::MolHDissociativeIzn;
  static constexpr auto type = "h2+ + e -> 2h+ + 2e";
};
RegisterComponent<H2PlusDissIzn> register_h2_plus_diss_izn;
struct D2PlusDissIzn : public hermes::MolHDissociativeIzn<hermes::HIsotope::d> {
  using MolHDissociativeIzn<hermes::HIsotope::d>::MolHDissociativeIzn;
  static constexpr auto type = "d2+ + e -> 2d+ + 2e";
};
RegisterComponent<D2PlusDissIzn> register_d2_plus_diss_izn;
struct T2PlusDissIzn : public hermes::MolHDissociativeIzn<hermes::HIsotope::t> {
  using MolHDissociativeIzn<hermes::HIsotope::t>::MolHDissociativeIzn;
  static constexpr auto type = "t2+ + e -> 2t+ + 2e";
};
RegisterComponent<T2PlusDissIzn> register_t2_plus_diss_izn;

// Elastic collisions between H isotope molecules and H isotope ions
struct H2HpElasticCollision : public hermes::MolHElasticCollision<hermes::HIsotope::h> {
  using MolHElasticCollision<hermes::HIsotope::h>::MolHElasticCollision;
  static constexpr auto type = "h2 + h+ -> h2 + h+";
};
RegisterComponent<H2HpElasticCollision> register_h2_hp_elastic_collision;
struct D2DpElasticCollision : public hermes::MolHElasticCollision<hermes::HIsotope::d> {
  using MolHElasticCollision<hermes::HIsotope::d>::MolHElasticCollision;
  static constexpr auto type = "d2 + d+ -> d2 + d+";
};
RegisterComponent<D2DpElasticCollision> register_d2_dp_elastic_collision;
struct T2TpElasticCollision : public hermes::MolHElasticCollision<hermes::HIsotope::t> {
  using MolHElasticCollision<hermes::HIsotope::t>::MolHElasticCollision;
  static constexpr auto type = "t2 + t+ -> t2 + t+";
};
RegisterComponent<T2TpElasticCollision> register_t2_tp_elastic_collision;

} // namespace

#endif // HYDROGEN_MOLECULE_REACTIONS_H
