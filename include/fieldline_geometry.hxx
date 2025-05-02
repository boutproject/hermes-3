#pragma once
#ifndef fieldline_geometry_H
#define fieldline_geometry_H

#include "component.hxx"
#include <bout/constants.hxx>

/**
 * @brief The ``fieldline_geometry`` component in Hermes-3 helps to set up 1D simulations and to compute cell geometry for other modules.
 */
struct FieldlineGeometry : public Component {

  /**
   * @brief Constructor for the FieldlineGeometry component.
   * @param name The name of the component.
   * @param options Options for the component.
   * @param solver Pointer to the solver.
   */
  FieldlineGeometry(std::string, Options& options, Solver*);

  /**
   * @brief Transform function (not implemented).
   * @param state Options containing the state.
   */
  void transform(Options& state) override;
  
  /**
   * @brief Output variables to the state.
   * @param state Options containing the state.
   */
  void outputVars(Options& state) override;

  private:
    bool diagnose;  ///< Flag to output additional diagnostics.

    Field3D lpar{0.0}; ///< Parallel distance from upstream.

    Field3D lambda_int{0.0}; ///< Radial width of the flux tube mapped upstream (lambda_q + 1.64S).
    Field3D fieldline_radius{0.0}; ///< Major radius of the magnetic fieldline.
    Field3D poloidal_magnetic_field{0.0}; ///< Poloidal magnetic field strength.
    Field3D toroidal_magnetic_field{0.0}; ///< Toroidal magnetic field strength along the fieldline.
    Field3D total_magnetic_field{0.0};   ///< Total magnetic field strength.
    Field3D pitch_angle{0.0};        ///< Pitch angle of the magnetic field (Bpol/B).
    
    Field3D transport_broadening{0.0}; ///< Flux tube broadening due to cross-field transport.
    Field3D flux_expansion{0.0};      ///< Flux expansion.

    Field3D flux_tube_width{0.0};    ///< Flux tube radial width.
    Field3D cell_poloidal_length{0.0}; ///< Poloidal length of the flux tube.
    Field3D cell_side_area{0.0};     ///< Poloidal length of the flux tube times its circumference.
    Field3D cell_volume{0.0};        ///< Poloidal length of the flux tube times its circumference and radial width.
};

namespace {
  RegisterComponent<FieldlineGeometry> register_fieldline_geometry("fieldline_geometry");
}

#endif // fieldline_geometry_H