#pragma once
#ifndef BRAGINSKII_CONDUCTION_H
#define BRAGINSKII_CONDUCTION_H

#include "component.hxx"

/// Calculates parallel heat conduction due to collisions
/// 
/// NOTE: This is global as that is the only way to ensure it gets run
/// after all collisions have been calculated. Logically this should
/// really apply species-by-species, but we can't do that until we
/// have dynamic ordering to make sure collision rates get calculated
/// first.
struct BraginskiiConduction : public Component {
  ///
  /// # Inputs
  ///
  /// - <component name>
  ///   - diagnose                     Output additional diagnostic fields?
  ///   - kappa_coefficient            Heat conduction constant. Default is 3.16 for
  ///                                  electrons, 3.9 otherwise
  ///   - kappa_limit_alpha            Flux limiter, off by default.
  ///   - conduction_collisions_mode   Can be multispecies: all collisions, or braginskii:
  ///                                  self collisions and ie
  ///
  /// - <species name>
  ///   - type                  Checks whether energy or pressure are evolved
  ///   - thermal_conduction    Include parallel heat conduction? Default is true
  ///
  BraginskiiConduction(std::string name, Options& alloptions, Solver*);


  /// Calculate conduction of energy for each species where this has been turned on.
  ///
  /// Uses
  ///   - species
  ///     - <name>
  ///       - AA
  ///       - collision_frequencies
  ///       - density
  ///       - temperature
  ///       - pressure
  ///
  /// Modifies
  ///   - species
  ///     - <name>
  ///       - energy_source     Conduction contribution to energy evolution
  ///       - kappa_par         The parallel heat conduction coefficient
  ///       - energy_flow_ylow  Energy flow diagnostics.
  ///
  void transform(Options &state) override;

  /// Add extra fields for output, or set attributes e.g docstrings
  void outputVars(Options &state) override;

private:
  std::map<std::string, Field3D> all_nu;   ///< Collision frequency for conduction
  std::map<std::string, Field3D> all_kappa_par; ///< Parallel heat conduction coefficient
  std::map<std::string, std::string> all_conduction_collisions_mode; ///< Collision selection, either multispecies or braginskii
  std::map<std::string, std::vector<std::string>> all_collision_names; ///< Collisions used for collisionality
  std::map<std::string, BoutReal> all_kappa_coefficient; ///< Leading numerical coefficient in parallel heat flux calculation
  std::map<std::string, BoutReal> all_kappa_limit_alpha; ///< Flux limit if >0
  std::map<std::string, Field3D> all_flow_ylow_conduction; ///< Conduction energy flow diagnostics
  /// Save more diagnostics?
  std::map<std::string, bool> all_diagnose;
};

namespace {
RegisterComponent<BraginskiiConduction> registercomponentbraginskiiconduction("braginskii_conduction");
}

#endif // BRAGINSKII_CONDUCTION_H
