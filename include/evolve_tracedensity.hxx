#pragma once
#ifndef EVOLVE_TRACEDENSITY_H
#define EVOLVE_TRACEDENSITY_H

#include "component.hxx"


struct EvolveTraceDensity : public Component {

  EvolveTraceDensity(std::string name, Options &options, Solver *solver);

  void transform(Options &state) override;

  void finally(const Options &state) override;

  void outputVars(Options &state) override;

  private:
  std::string name;     ///< Short name of species e.g "e"

  BoutReal charge;      ///< Species charge e.g. electron = -1
  BoutReal AA;          ///< Atomic mass e.g. proton = 1
  
  Field3D tN;            ///< Species density (normalised, evolving)

  bool bndry_flux;      ///< Allow flows through boundaries?
  bool exb_advection;   ///< Include ExB advection?
  bool poloidal_flows;  ///< Include ExB flow in Y direction?
  bool neumann_boundary_average_z; ///< Apply neumann boundary with Z average?
  bool disable_ddt;
  bool dissipative;
  BoutReal density_floor;
  bool low_n_diffuse;   ///< Parallel diffusion at low density
  bool low_n_diffuse_perp;  ///< Perpendicular diffusion at low density
  BoutReal pressure_floor; ///< When non-zero pressure is needed
  bool low_p_diffuse_perp; ///< Add artificial cross-field diffusion at low pressure?
  BoutReal hyper_z;    ///< Hyper-diffusion in Z
  BoutReal hyper_n;
  BoutReal scale_ExB;
  bool evolve_log; ///< Evolve logarithm of density?
  Field3D logN;    ///< Logarithm of density (if evolving)
  bool isMMS;
  Field3D source, final_source; ///< External input source
  Field3D Sn; ///< Total density source
  bool sourceterm;
  BoutReal adapt_source;
  
  bool source_only_in_core;  ///< Zero source where Y is non-periodic?
  bool source_time_dependent; ///< Is the input source time dependent?
  BoutReal source_normalisation; ///< Normalisation factor [m^-3/s]
  BoutReal time_normalisation; ///< Normalisation factor [s]
  FieldGeneratorPtr source_prefactor_function;

  /// Modifies the `source` member variable
  void updateSource(BoutReal time);
  bool output_ddt;
  bool diagnose; ///< Output additional diagnostics?
  Field3D flow_xlow, flow_ylow; ///< Particle flow diagnostics

  Coordinates::FieldMetric bracket_factor; ///< For non-Clebsch coordinate systems (e.g. FCI)
};

namespace {
RegisterComponent<EvolveTraceDensity> registercomponentevolvedensity("evolve_tracedensity");
}


#endif // EVOLVE_TRACEDENSITY_H
