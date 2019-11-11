#pragma once

#ifndef __ATOMICPP_PRAD_H__
#define __ATOMICPP_PRAD_H__

#include "ImpuritySpecies.hxx"

/// Calculates the relative distribution across ionisation stages of the
/// impurity by assuming collisional-radiative
/// equilbrium. This is then used to calculate the density within each state,
/// allowing the total power at a point
/// to be evaluated
///
/// Inputs
/// ------
///
/// @param[in] impurity   An object describing a species
/// @param[in] Te         Electron temperature in eV
/// @param[in] Ne         Electron density in m^-3
/// @param[in] Ni         Impurity density in m^-3
///                       summed over all charge states
/// @param[in] Nn         Neutral atomic density in m^-3
///
/// Returns
/// -------
/// Radiated power in W/m^-3. Since this routine assumes
/// collisional-radiative equilibrium, then this power is the
/// same as the electron cooling rate.
/// 
/// Example
/// -------
///
/// ImpuritySpecies impurity("c"); // Carbon
///
/// double total_power = computeRadiatedPower(impurity, Te, Ne, Ni, Nn);
///
double computeRadiatedPower(ImpuritySpecies &impurity, double Te, double Ne,
                            double Ni, double Nn);

#endif // __ATOMICPP_PRAD_H__

