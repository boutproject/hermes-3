#pragma once
#ifndef BRAGINSKII_CONDUCTION_H
#define BRAGINSKII_CONDUCTION_H

#include <map>
#include <string>
#include <vector>

#include <bout/bout_types.hxx>
#include <bout/field3d.hxx>
#include <bout/options.hxx>

#include "component.hxx"

/// Calculates parallel heat conduction due to collisions
///
/// NOTE: This is global as that is the only way to ensure it gets run
/// after all collisions have been calculated. Logically this should
/// really apply species-by-species, but we can't do that until we
/// have dynamic ordering to make sure collision rates get calculated
/// first.
struct BraginskiiConduction : public NamedComponent<BraginskiiConduction> {
  ///
  /// # Inputs
  ///
  /// - <component name>
  ///   - diagnose                     Output additional diagnostic fields?
  ///   - kappa_coefficient            Heat conduction constant. Default is 3.16 for
  ///                                  electrons, 3.9 otherwise
  ///   - kappa_limit_alpha            Flux limiter, off by default. Free-streaming coefficient.
  ///   - kappa_limit_model            Limiter model: "local" or "connection_length".
  ///   - kappa_limit_q95              Safety factor used by the connection-length limiter.
  ///   - kappa_limit_R                Major radius used by the connection-length limiter [m].  
  /// 
  ///   - conduction_collisions_mode   Can be multispecies: all collisions, or braginskii:
  ///                                  self collisions and ie
  ///
  /// - <species name>
  ///   - type                  Checks whether energy or pressure are evolved
  ///   - thermal_conduction    Include parallel heat conduction? Default is true
  ///
  BraginskiiConduction(const std::string& name, Options& alloptions, Solver*);

  /// Add extra fields for output, or set attributes e.g docstrings
  void outputVars(Options& state) override;

  static constexpr auto type = "braginskii_conduction";

private:
  std::map<std::string, Field3D> all_nu;        ///< Collision frequency for conduction
  std::map<std::string, Field3D> all_kappa_par; ///< Parallel heat conduction coefficient
  std::map<std::string, std::string>
      all_conduction_collisions_mode; ///< Collision selection, either multispecies or
                                      ///< braginskii
  std::map<std::string, std::vector<std::string>>
      all_collision_names; ///< Collisions used for collisionality
  std::map<std::string, BoutReal>
      all_kappa_coefficient; ///< Leading numerical coefficient in parallel heat flux
                             ///< calculation
  std::map<std::string, BoutReal> all_kappa_limit_alpha; ///< Flux limit if >0

  std::map<std::string, std::string>
      all_kappa_limit_model; ///< "local" or "connection_length"

  std::map<std::string, BoutReal>
      all_kappa_limit_q95; ///< Safety factor used by connection-length limiter

  std::map<std::string, BoutReal>
      all_kappa_limit_R; ///< Major radius, normalised to the Hermes length unit rho_s0

  std::map<std::string, Field3D>
      all_flow_ylow_conduction; ///< Conduction energy flow diagnostics
  /// Save more diagnostics?
  std::map<std::string, bool> all_diagnose;

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
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<BraginskiiConduction> registercomponentbraginskiiconduction;
}

#endif // BRAGINSKII_CONDUCTION_H
