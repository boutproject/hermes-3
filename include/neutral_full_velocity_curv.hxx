
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
  bool disable_dndt;
  bool include_D, include_nu;
  Field3D anomalous_D, anomalous_nu;
  bool inherited_T;
  bool parallel_dirichlet;
  bool neutral_viscosity; ///< include viscosity?
  bool neutral_conduction; ///< Include heat conduction?
  bool evolve_momentum; ///< Evolve parallel momentum?
  bool evolve_momentum_xz;
  bool evolve_pressure;
  bool momentum_advection;
  bool momentum_loss;
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


  struct Stencil1D {
    // Cell centre values                                                                                                                            
    BoutReal c, m, p, mm, pp;
    // Left and right cell face values                                                                                                              
    BoutReal L, R;
  };

  BoutReal minmod(BoutReal a, BoutReal b) {
    if (a * b <= 0.0)
      return 0.0; 
    if (fabs(a) < fabs(b))
      return a;
    return b;
  }

  BoutReal minmod(BoutReal a, BoutReal b, BoutReal c) {
    // If any of the signs are different, return zero gradient                                                                                      
    if ((a * b <= 0.0) || (a * c <= 0.0)) {
      return 0.0;
    }
    // Return the minimum absolute value                                                                                                            
    return SIGN(a) * BOUTMIN(fabs(a), fabs(b), fabs(c));
  }

  void MinMod(Stencil1D& n, const BoutReal h) {
    // Choose the gradient within the cell                                                                                                           
    // as the minimum (smoothest) solution                                                                                                           
    BoutReal slope = minmod(n.p - n.c, n.c - n.m);
    n.L = n.c - 0.5 * slope; // 0.25*(n.p - n.m);
    n.R = n.c + 0.5 * slope; // 0.25*(n.p - n.m)
  }

  // Monotonized Central limiter (Van-Leer)                                                                                                             
  void MC(Stencil1D& n, const BoutReal h) {
    BoutReal slope = minmod(2. * (n.p - n.c), 0.5 * (n.p - n.m), 2. * (n.c - n.m));
    n.L = n.c - 0.5 * slope;
    n.R = n.c + 0.5 * slope;
  }



  const Field3D  Div_perp(Field3D& f, Field3D& vx, Field3D& vz, Field3D& spd, const bool& dissipation=false) {

    Mesh* mesh = f.getMesh();
    Coordinates* coord = mesh->getCoordinates();

    Field3D cellcross_x = coord->cellvolume / (coord->dx * sqrt(coord->g_11));
    Field3D cellcross_z = coord->cellvolume / (coord->dz * sqrt(coord->g_33));

    Field3D result{zeroFrom(f)};

    BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
      const auto ixp = i.xp();
      const auto ixpp = i.xpp();
      const auto ixm = i.xm();
      const auto ixmm = i.xmm();
      
      const auto izp = i.zp();
      const auto izpp = i.zpp();
      const auto izm = i.zm();
      const auto izmm = i.zmm();

      BoutReal cellarea_xup = 0.5 * (cellcross_x[i] + cellcross_x[ixp]);
      BoutReal cellarea_xdown = 0.5 * (cellcross_x[i] + cellcross_x[ixm]);

      BoutReal cellarea_zup = 0.5 * (cellcross_z[i] + cellcross_z[izp]);
      BoutReal cellarea_zdown = 0.5 * (cellcross_z[i] + cellcross_z[izm]);

      BoutReal v_xup = 0.5 * (vx[i] + vx[ixp]);
      BoutReal v_xdown = 0.5 * (vx[i] + vx[ixm]);

      BoutReal v_zup = 0.5 * (vz[i] + vz[izp]);
      BoutReal v_zdown = 0.5 * (vz[i] + vz[izm]);

      ////////////////////////////////////////////////////////////////////////7
      //                         x- direction
      
      Stencil1D sx;
      sx.c = f[i];
      sx.m = f[ixm];
      sx.p = f[ixp];      
      MinMod(sx, coord->dx[i]);

      Stencil1D sxm;
      sxm.c = f[ixm];
      sxm.m = f[ixmm];
      sxm.p = f[i];
      MinMod(sxm, coord->dx[i]);

      Stencil1D sxp;
      sxp.c = f[ixp];
      sxp.m = f[i];
      sxp.p = f[ixpp];
      MinMod(sxp, coord->dx[i]);
      
      BoutReal amax_xup = BOUTMAX(fabs(vx[i]),fabs(vx[ixp]),spd[i],spd[ixp]);
      BoutReal flux_xup = (0.5 * (sx.R + sxp.L) * v_xup + 0.5 * amax_xup * (sx.R - sxp.L)) * cellarea_xup;

      BoutReal amax_xdown = BOUTMAX(fabs(vx[i]),fabs(vx[ixm]),spd[i],spd[ixm]);
      BoutReal flux_xdown = (0.5 * (sx.L + sxm.R) * v_xdown + 0.5 * amax_xdown * (sxm.R - sx.L)) * cellarea_xdown;
      
      

      /////////////////////////////////////////////////////////
      //         Z direction

      Stencil1D sz;
      sz.c = f[i];
      sz.m = f[izm];
      sz.p = f[izp];
      MinMod(sz, coord->dz[i]);

      Stencil1D szm;
      szm.c = f[izm];
      szm.m = f[izmm];
      szm.p = f[i];
      MinMod(szm, coord->dz[i]);

      Stencil1D szp;
      szp.c = f[izp];
      szp.m = f[i];
      szp.p = f[izpp];
      MinMod(szp, coord->dz[i]);
      
      
      BoutReal amax_zup = BOUTMAX(fabs(vz[i]),fabs(vz[izp]),spd[i],spd[izp]);
      BoutReal flux_zup = (0.5 * (sz.R + szp.L) * v_zup + 0.5 * amax_zup * (sz.R - szp.L)) * cellarea_zup;

      BoutReal amax_zdown = BOUTMAX(fabs(vz[i]),fabs(vz[izm]),spd[i],spd[izm]);
      BoutReal flux_zdown = (0.5 * (sz.L + szm.R) * v_zdown + 0.5 * amax_zdown * (szm.R - sz.L)) * cellarea_zdown;

      result[i] = (flux_xup + flux_zup - flux_xdown - flux_zdown) / coord->cellvolume[i];
      
    }

    return result;
  }

  
  
  const Field3D Div_perp_flows(Field3D& f, Field3D& v, Field3D& vx, Field3D& vz, Field3D& spd) {
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

      BoutReal vx_up = 0.5 * (vx[i] + vx[ixp]);
      BoutReal vx_down = 0.5 * (vx[i] + vx[ixm]);

      BoutReal vz_up = 0.5 * (vz[i] + vz[izp]);
      BoutReal vz_down = 0.5 * (vz[i] + vz[izm]);

      BoutReal v_xup = (0.5 * (v[i] + v[ixp]));
      BoutReal v_xdown = (0.5 * (v[i] + v[ixm]));
      BoutReal v_zup = (0.5 * (v[i] + v[izp]));
      BoutReal v_zdown = (0.5 * (v[i] + v[izm]));
      
      BoutReal flux_xup = f_xup * vx_up * v_xup * cellarea_xup;
      BoutReal flux_xdown = f_xdown * vx_down * v_xdown * cellarea_xdown;

      BoutReal flux_zup = f_zup * vz_up * v_zup * cellarea_zup;
      BoutReal flux_zdown = f_zdown * vz_down * v_zdown * cellarea_zdown;

      result[i] = (flux_xup + flux_zup - flux_xdown -flux_zdown) / coord->cellvolume[i];
    }

    return result;

    
  }

  
  Field3D Div_a_Grad_perp(Field3D a, Field3D b, bool use_finite=false) {
    if (a.isFci()) {
      if (use_finite) {
	return Div_a_Grad_perp_curv(a,b);
      } else {
	return (*dagp)(a, b, false);
      }
    }
    return FV::Div_a_Grad_perp(a, b);
  }
};

namespace {
RegisterComponent<NeutralFullVelocityCurv> registersolverneutralmixed("neutral_full_velocity_curv");
}

#endif // NEUTRAL_FULL_VELOCITY_CURV_H
