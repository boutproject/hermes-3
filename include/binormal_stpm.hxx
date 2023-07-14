#pragma once
#ifndef BINORMAL_STPM_H
#define BINORMAL_STPM_H

#include <bout/field3d.hxx>

#include "component.hxx"

/// Adds terms to Density, Pressure and Momentum equations following the
/// stellarator 2-point model from Yuhe Feng et al., PPCF 53 (2011) 024009
/// The terms add the effective parallel contribution of perpendicular transport
/// in long connection length scenarios.
/// B Shanahan 2023 <brendan.shanahan@ipp.mpg.de>
struct BinormalSTPM : public Component {
  ///
  /// # Inputs
  ///
  /// - <name>
  ///   - D           Perpendicular density diffusion coefficient
  ///   - chi         Perpendicular heat diffusion coefficient
  ///   - nu          Perpendicular momentum diffusion coefficient
  ///   - Vbn         Binormal velocity from cross-field drifts
  ///   - Theta       Field line pitch as described by Feng et al.
  ///
  BinormalSTPM(std::string name, Options& options, Solver* solver);
  
  /// Sets
  /// - species
  ///   - <name>
  ///     - pressure correction 
  ///     - momentum correction
  ///     - density correction
  ///
  void transform(Options& state) override;


private:
  std::string name; ///< Short name of the species e.g. h+
  Field3D Theta, chi, D, nu; ///< Field line pitch, anomalous thermal, momentum diffusion
  Field3D Vbn;  ///< Binormal velocity from cross-field drifts input
  Field3D nu_Theta, chi_Theta, D_Theta; ///< nu/Theta, chi/Theta, D/Theta, precalculated

};

namespace {
RegisterComponent<BinormalSTPM> registercomponentbinormalstpm("binormal_stpm");
}

#endif // BINORMAL_STPM
