#pragma once
#ifndef BRAGINSKII_THERMAL_FORCE_H
#define BRAGINSKII_THERMAL_FORCE_H

#include <string>

#include <bout/options.hxx>

#include "component.hxx"

/// Simple calculation of the thermal force for the Braginskii closure
///
/// Important: This implements a quite crude approximation,
/// which is intended for initial development and testing.
/// The expressions used are only valid for trace heavy ions and
/// light main ion species, and would not be valid for Helium impurities
/// in a D-T plasma, for example. For this reason only collisions
/// where one ion has an atomic mass < 4, and the other an atomic mass > 10
/// are considered. Warning messages will be logged for species combinations
/// which are not calculated.
///
/// Options used:
///
/// - <name>
///   - electron_ion  : bool   Include electron-ion collisions?
///   - ion_ion       : bool   Include ion-ion elastic collisions?
///
struct BraginskiiThermalForce : public Component {
  BraginskiiThermalForce(const std::string& name, Options& alloptions, Solver*) {
    Options& options = alloptions[name];
    electron_ion = options["electron_ion"]
                       .doc("Include electron-ion collisions?")
                       .withDefault<bool>(true);

    ion_ion = options["ion_ion"]
                  .doc("Include ion-ion elastic collisions?")
                  .withDefault<bool>(true);
  }

private:
  bool electron_ion; ///< Include electron-ion collisions?
  bool ion_ion; ///< Include ion-ion elastic collisions?

  bool first_time{true}; ///< True the first time transform() is called

  /// Inputs
  /// - species
  ///   - e           [ if electron_ion true ]
  ///     - charge
  ///     - density
  ///     - temperature
  ///   - <species>
  ///     - charge    [ Checks, skips species if not set ]
  ///     - AA
  ///     - temperature [ If AA < 4 i.e. "light" species ]
  ///
  /// Outputs
  /// - species
  ///   - e
  ///     - momentum_source  [ if electron_ion true ]
  ///   - <species>          [ if AA < 4 ("light") or AA > 10 ("heavy") ]
  ///     - momentum_source
  ///
  void transform(GuardedOptions &state) override;
};

namespace {
RegisterComponent<BraginskiiThermalForce>
    registercomponentbraginskiithermalforce("thermal_force");
}

#endif // BRAGINSKII_THERMAL_FORCE_H
