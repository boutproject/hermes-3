
#pragma once
#ifndef NEUTRAL_FULL_VELOCITY_CURV_H
#define NEUTRAL_FULL_VELOCITY_CURV_H

#include <memory>
#include <string>

#include <bout/invert_laplace.hxx>

#include "component.hxx"
#include <bout/yboundary_regions.hxx>



/// Evolve density, parallel momentum and pressure
/// for a neutral gas species with cross-field diffusion
struct NeutralFullVelocityCurv : public Component {
  ///
  /// @param name     The name of the species e.g. "h"
  /// @param options  Top-level options. Settings will be taken from options[name]
  /// @param solver   Time-integration solver to be used
  NeutralFullVelocityCurv(const std::string& name, Options& options, Solver *solver);
  
  /// Modify the given simulation state
  void transform(Options &state) override;
  
  /// Use the final simulation state to update internal state
  /// (e.g. time derivatives)
  void finally(const Options &state) override;

  /// Add extra fields for output, or set attributes e.g docstrings
  void outputVars(Options &state) override;

  /// Preconditioner
  void precon(const Options &state, BoutReal gamma) override;
private:
  std::string name;  ///< Species name
  YBoundary yboundary;
  std::shared_ptr<FCI::dagp_fv> dagp;
  Field3D Nn, Pn, NVn; // Density, pressure and parallel momentum
  Field3D NVn_x, NVn_z;
  Field3D Vn, Vn_x, Vn_z; ///< Neutral parallel velocity
  Field3D Tn; ///< Neutral temperature
  Field3D Nnlim, Pnlim, logPnlim, Vnlim, Tnlim; // Limited in regions of low density
  bool isMMS;
  BoutReal AA; ///< Atomic mass (proton = 1)
  BoutReal n_lowsource, T_lowsource, lowsource_scale;

  Field3D Rnn;
  Field3D Dnn;
  BoutReal neutral_lmax;
  
  BoutReal temperature_floor;

  bool dissipative;
  BoutReal density_floor; ///< Minimum Nn used when dividing NVn by Nn to get Vn.
  BoutReal pressure_floor; ///< Minimum Pn used when dividing Pn by Nn to get Tn.


  bool parallel_dirichlet;
  bool neutral_viscosity; ///< include viscosity?
  bool neutral_conduction; ///< Include heat conduction?
  bool evolve_momentum; ///< Evolve parallel momentum?
  bool evolve_momentum_xz;
  bool evolve_pressure;
  Field3D initial_Tn;
  Field3D initial_Vn, initial_Vn_x,initial_Vn_z;

  
  bool use_finite_difference;
  Field3D kappa_n, eta_n; ///< Neutral conduction and viscosity

  bool precondition {true}; ///< Enable preconditioner?
  bool lax_flux; ///< Use Lax flux for advection terms
  std::unique_ptr<Laplacian> inv; ///< Laplacian inversion used for preconditioning

  Field3D density_source, pressure_source; ///< External input source
  Field3D Sn, Sp, Snv; ///< Particle, pressure and momentum source
  Field3D sound_speed; ///< Sound speed for use with Lax flux

  bool output_ddt; ///< Save time derivatives?
  bool diagnose; ///< Save additional diagnostics?

  
  // Flow diagnostics
  Field3D pf_adv_perp_xlow, pf_adv_perp_ylow, pf_adv_par_ylow;
  Field3D mf_adv_perp_xlow, mf_adv_perp_ylow, mf_adv_par_ylow;
  Field3D mf_visc_perp_xlow, mf_visc_perp_ylow, mf_visc_par_ylow;
  Field3D ef_adv_perp_xlow, ef_adv_perp_ylow, ef_adv_par_ylow;
  Field3D ef_cond_perp_xlow, ef_cond_perp_ylow, ef_cond_par_ylow;



  const Field3D Grad_x(Field3D& a) {
    Mesh* mesh = a.getMesh();
    Coordinates* coord = mesh->getCoordinates();
    return DDX(a) / sqrt(coord->g_11);
  }

  
  const Field3D Grad_z(Field3D& a) {
    Mesh* mesh = a.getMesh();
    Coordinates* coord = mesh->getCoordinates();
    return DDZ(a) / sqrt(coord->g_33);
  }


  const Field3D  Div_perp(Field3D& f, Field3D& vx, Field3D& vz, Field3D& spd) {

    Mesh* mesh = f.getMesh();
    Coordinates* coord = mesh->getCoordinates();

    Field3D cellcross_x = coord->cellvolume / (coord->dx * sqrt(coord->g_11));
    Field3D cellcross_z = coord->cellvolume / (coord->dz * sqrt(coord->g_33));

    Field3D result{zeroFrom(f)};

    BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
      const auto ixp = i.xp();
      const auto ixm = i.xm();

      const auto izp = i.zp();
      const auto izm = i.zm();

      BoutReal cellarea_xup = 0.5 * (cellcross_x[i] + cellcross_x[ixp]);
      BoutReal cellarea_xdown = 0.5 * (cellcross_x[i] + cellcross_x[ixm]);

      BoutReal cellarea_zup = 0.5 * (cellcross_z[i] + cellcross_z[izp]);
      BoutReal cellarea_zdown = 0.5 * (cellcross_z[i] + cellcross_z[izm]);
      
      BoutReal f_xup = 0.5 * (f[i] + f[ixp]);
      BoutReal f_xdown = 0.5 * (f[i] + f[ixm]);

      BoutReal f_zup = 0.5 * (f[i] + f[izp]);
      BoutReal f_zdown = 0.5 * (f[i] + f[izm]);

      BoutReal v_xup = 0.5 * (vx[i] + vx[ixp]);
      BoutReal v_xdown = 0.5 * (vx[i] + vx[ixm]);

      BoutReal v_zup = 0.5 * (vz[i] + vz[izp]);
      BoutReal v_zdown = 0.5 * (vz[i] + vz[izm]);

      BoutReal flux_xup = f_xup * v_xup * cellarea_xup;
      BoutReal flux_xdown = f_xdown * v_xdown * cellarea_xdown;

      BoutReal flux_zup	= f_zup * v_zup * cellarea_zup;
      BoutReal flux_zdown = f_zdown * v_zdown * cellarea_zdown;

      result[i] = (flux_xup + flux_zup - flux_xdown -flux_zdown) / coord->cellvolume[i];            
    }

    return result;
  }

  
  const Field3D Div_perp_fvv_x(Field3D& f,Field3D& v, Field3D& spd) {
    Field3D result{zeroFrom(f)};
  }

  const Field3D Div_perp_fvv_z(Field3D& f,Field3D& v, Field3D& spd) {
    Field3D result{zeroFrom(f)};
  }
  


  
  Field3D Div_a_Grad_perp(Field3D a, Field3D b) {
    if (a.isFci()) {
      return (*dagp)(a, b, false);
    }
    return FV::Div_a_Grad_perp(a, b);
  }
};

namespace {
RegisterComponent<NeutralFullVelocityCurv> registersolverneutralmixed("neutral_full_velocity_curv");
}

#endif // NEUTRAL_FULL_VELOCITY_CURV_H
