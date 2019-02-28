/*
 * Base class for neutral gas model
 */

#ifndef __NEUTRAL_MODEL_H__
#define __NEUTRAL_MODEL_H__

#include <field3d.hxx>
#include <options.hxx>
#include <bout/solver.hxx>
#include <bout/mesh.hxx>

#include "radiation.hxx"

class NeutralModel {
public:
  NeutralModel(Options &options) {
    S = 0;
    F = 0;
    Fperp = 0;
    Qi = 0;
    Rp = 0;
    Rn = 0;

    SAVE_REPEAT4(S,F,Qi,Rp);
    
    // Options for calculating rates
    OPTION(options, Eionize, 30); // Energy loss per ionisation [eV]
  }
  virtual ~NeutralModel() {}
  
  /*!
   * Creates an instance of NeutralModel, based on given options
   */
  static NeutralModel* create(Solver *solver, Mesh *mesh, Options &options);
  
  /*!
   * Set normalisations for temperature [eV], density [m^-3], 
   * length [m] and frequency [s^-1]
   * 
   */
  void setNormalisation(BoutReal Te, BoutReal Ne, BoutReal B, BoutReal length, BoutReal freq) { 
    Tnorm = Te; Nnorm = Ne; Bnorm = B; Lnorm = length; Fnorm = freq;
  }

  /*!
   * Update plasma quantities
   */
  virtual void update(const Field3D &Ne, const Field3D &Te, const Field3D &Ti, const Field3D &Vi) = 0;
  
  /*!
   *
   */
  virtual void addDensity(int UNUSED(x), int UNUSED(y), int UNUSED(z), BoutReal UNUSED(dndt)) {}
  virtual void addPressure(int UNUSED(x), int UNUSED(y), int UNUSED(z), BoutReal UNUSED(dpdt)) {}
  virtual void addMomentum(int UNUSED(x), int UNUSED(y), int UNUSED(z), BoutReal UNUSED(dnvdt)) {}
  
  /*!
   * Preconditioning
   */
  virtual void precon(BoutReal UNUSED(t), BoutReal UNUSED(gamma), BoutReal UNUSED(delta)) {}

  Field3D S;     // Plasma particle sink, neutral source
  Field3D Qi;    // Power transfer from ions
  Field3D F;     // Ion-neutral friction
  Field3D Fperp; // Ion-neutral friction in vorticity
  Field3D Rp;    // Radiation from the plasma
  Field3D Rn;    // Radiation from the neutrals

protected:
  BoutReal Tnorm, Nnorm, Bnorm, Lnorm, Fnorm; // Normalisations for temperature, density, magnetic field, lengths and frequencies
  
  UpdatedRadiatedPower hydrogen; // Atomic rates (H.Willett)
  
  BoutReal Eionize;   // Energy loss per ionisation [eV]
  
  void neutral_rates(const Field3D &Ne, const Field3D &Te, const Field3D &Ti, const Field3D &Vi,    // Plasma quantities
                     const Field3D &Nn, const Field3D &Tn, const Field3D &Vnpar, // Neutral gas
                     Field3D &S, Field3D &F, Field3D &Qi, Field3D &R,  // Transfer rates
                     Field3D &Riz, Field3D &Rrc, Field3D &Rcx);        // Rates

private:
  NeutralModel();
};

#endif // __NEUTRAL_MODEL_H__
