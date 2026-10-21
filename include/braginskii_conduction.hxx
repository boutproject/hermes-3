#pragma once
#ifndef BRAGINSKII_CONDUCTION_H
#define BRAGINSKII_CONDUCTION_H

#include <map>
#include <string>
#include <vector>

#include <bout/bout_types.hxx>
#include <bout/field3d.hxx>
#include <bout/options.hxx>

#if BOUT_HAS_PETSC
#include <bout/petsc_operators.hxx>
#endif

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
  ///   - kappa_limit_alpha            Flux limiter, off by default.
  ///   - kappa_free_streaming         Reduce kappa by 1 / (1 + lambda / (c L_par)),
  ///                                  lambda = v_th tau, 1/L_par = 1/L_line (+ |Grad_par T|/T).
  ///                                  Default false
  ///   - kappa_free_streaming_c       c, default 1
  ///   - kappa_free_streaming_gradient  Include |Grad_par T|/T in 1/L_par? Default false
  ///   - kappa_density_threshold      n_th [m^-3]: kappa -> kappa n / (n + n_cell), with
  ///                                  n_cell = n_th (ds_ref / ds)^2 for the parallel cell
  ///                                  length ds. At low density the cell conduction rate
  ///                                  chi / ds^2 is then bounded by the same value in every
  ///                                  cell; the limiter vanishes for n >> n_cell.
  ///                                  Default -1 (off)
  ///   - kappa_density_threshold_length  ds_ref [m]. Default: mean cell length along
  ///                                  each field line (needs a grid without branch cuts)
  ///   - kappa_cell_timescale         Minimum conduction time [s] across one parallel
  ///                                  cell: kappa -> kappa / (1 + kappa / kappa_max),
  ///                                  kappa_max = 1.5 n ds^2 / timescale. Default -1 (off)
  ///   - kappa_free_streaming_length  L_line [m]. Default: field line length from the
  ///                                  grid (needs a grid without branch cuts)
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
  /// Free-streaming reduction of kappa (species with kappa_free_streaming = true):
  /// coefficient c and field line length
  std::map<std::string, BoutReal> all_kappa_free_streaming_c;
  std::map<std::string, Field2D> all_L_line;
  std::map<std::string, bool> all_kappa_free_streaming_gradient;
  /// Low-density limiter threshold per cell, n_th (ds_ref / ds)^2 (normalised);
  /// absent if off
  std::map<std::string, Field2D> all_kappa_density_threshold;
  /// kappa_max / n = 1.5 ds^2 / timescale of the cell timescale cap; absent if off
  std::map<std::string, Field2D> all_kappa_cell_max_per_n;
  std::map<std::string, Field3D>
      all_flow_ylow_conduction; ///< Conduction energy flow diagnostics
  /// Save more diagnostics?
  std::map<std::string, bool> all_diagnose;

  bool use_petsc_operator; ///< Use Petsc operators for parallel heat conduction?

#if BOUT_HAS_PETSC
  PetscOperators::Parallel
      petsc_parallel; ///< Petsc parallel operator object used for conduction
#endif

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

  bout::ConductionMethod conduction_method{bout::ConductionMethod::Original};
};

namespace {
RegisterComponent<BraginskiiConduction> registercomponentbraginskiiconduction;
}

#endif // BRAGINSKII_CONDUCTION_H
