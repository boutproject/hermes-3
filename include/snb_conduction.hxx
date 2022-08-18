#pragma once
#ifndef SNB_CONDUCTION_H
#define SNB_CONDUCTION_H

#include "component.hxx"

#include <bout/snb.hxx>

/// Calculate electron heat flux using the Shurtz-Nicolai-Busquet (SNB) model
///
/// This component will only calculate divergence of heat flux for the
/// electron (`e`) species.
///
/// # Usage
///
/// Add as a top-level component after both electron temperature and
/// collision times have been calculated.
///
/// Important: If evolving electron pressure, disable thermal
/// conduction or that will continue to add Spitzer heat conduction.
///
/// ```
/// [hermes]
/// components = e, ..., collisions, snb_conduction
///
/// [e]
/// type = evolve_pressure, ...
/// thermal_conduction = false # For evolve_pressure
///
/// [snb_conduction]
/// diagnose = true # Saves heat flux diagnostics
/// ```
///
/// # Useful references:
///
///  *  Braginskii equations by R.Fitzpatrick:
///     http://farside.ph.utexas.edu/teaching/plasma/Plasmahtml/node35.html
///
///  *  J.P.Brodrick et al 2017: https://doi.org/10.1063/1.5001079 and
///     https://arxiv.org/abs/1704.08963
///
///  *  Shurtz, Nicolai and Busquet 2000: https://doi.org/10.1063/1.1289512
///
struct SNBConduction : public Component {

  /// Inputs
  ///  - <name>
  ///    - diagnose   Saves Div_Q_SH and Div_Q_SNB
  SNBConduction(std::string name, Options& alloptions, Solver*) : snb(alloptions[name]) {
    AUTO_TRACE();
    auto& options = alloptions[name];

    if (options["diagnose"]
        .doc("Save additional output diagnostics")
          .withDefault<bool>(false)) {
      SAVE_REPEAT(Div_Q_SH, Div_Q_SNB);
    }
  }

  /// Inputs
  /// - species
  ///   - e
  ///     - density
  ///     - collision_frequency
  ///
  /// Sets
  /// - species
  ///   - e
  ///     - energy_source
  void transform(Options& state) override;

private:
  bout::HeatFluxSNB snb;

  Field3D Div_Q_SH, Div_Q_SNB; ///< Divergence of heat fluxes
};

namespace {
RegisterComponent<SNBConduction> registercomponentsnbconduction("snb_conduction");
}

#endif // SNB_CONDUCTION_H

