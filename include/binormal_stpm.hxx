#pragma once
#ifndef BINORMAL_STPM_H
#define BINORMAL_STPM_H

#include <bout/field3d.hxx>

#include "component.hxx"

/// Adds terms to Density, Pressure and Momentum equations following the
/// stellarator 2-point model from Yuhe Feng et al., PPCF 53 (2011) 024009
/// The terms add the effective parallel contribution of perpendicular transport
/// which is of significance in long connection length scenarios.
/// B Shanahan 2023 <brendan.shanahan@ipp.mpg.de>
struct BinormalSTPM : public Component {
  ///
  /// # Inputs
  ///
  /// - <name>
  ///   - D           Perpendicular density diffusion coefficient
  ///   - chi         Perpendicular heat diffusion coefficient
  ///   - nu          Perpendicular momentum diffusion coefficient
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
  void outputVars(Options &state) override;


private:
  std::string name; ///< Short name of the species e.g. h+
  bool diagnose; ///< Output diagnostics?
  Field3D Theta, chi, D, nu; ///< Field line pitch, anomalous thermal, momentum diffusion
  Field3D nu_Theta, chi_Theta, D_Theta; ///< nu/Theta, chi/Theta, D/Theta, precalculated
  Field3D Theta_inv; ///< Precalculate 1/Theta

};

namespace {
RegisterComponent<BinormalSTPM> registercomponentbinormalstpm("binormal_stpm");
}

#endif // BINORMAL_STPM
