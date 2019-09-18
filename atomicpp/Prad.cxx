// Program name: atomic++/Prad.cpp
// Author: Thomas Body
// Author email: tajb500@york.ac.uk
// Date of creation: 17 July 2017
//
// Program function: output the radiated power (Prad)
//                   by using OpenADAS rates on output JSON from SD1D run
//
// Based on the TBody/atomic1D code, which is in turn based on the cfe316/atomic
// code
//

// Include declarations
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept> //For error-throwing
#include <string>
#include <vector>

#include <memory>

#include "Prad.hxx"
#include "RateCoefficient.hxx"

using namespace std;

double computeRadiatedPower(ImpuritySpecies &impurity, double Te, double Ne,
                            double Ni, double Nn) {
  
  int Z = impurity.get_atomic_number();
  vector<double> iz_stage_distribution(Z + 1, 0.0);

  // Set GS density equal to 1 (arbitrary)
  iz_stage_distribution[0] = 1;
  double sum_iz = 1;

  // Loop over 0, 1, ..., Z-1
  // Each charge state is set in terms of the density of the previous
  for (int k = 0; k < Z; ++k) {
    // Ionisation
    // Get the RateCoefficient from the rate_coefficient map (atrribute of
    // impurity)
    shared_ptr<RateCoefficient> iz_rate_coefficient =
        impurity.get_rate_coefficient("ionisation");
    // Evaluate the RateCoefficient at the point
    double k_iz_evaluated = iz_rate_coefficient->call0D(k, Te, Ne);

    // Recombination
    // Get the RateCoefficient from the rate_coefficient map (atrribute of
    // impurity)
    shared_ptr<RateCoefficient> rec_rate_coefficient =
        impurity.get_rate_coefficient("recombination");
    // Evaluate the RateCoefficient at the point
    double k_rec_evaluated = rec_rate_coefficient->call0D(k, Te, Ne);

    // The ratio of ionisation from the (k)th stage and recombination from the
    // (k+1)th sets the equilibrium densities
    // of the (k+1)th stage in terms of the (k)th (since R = Nz * Ne *
    // rate_coefficient) N.b. Since there is no
    // ionisation from the bare nucleus, and no recombination onto the neutral
    // (ignoring anion formation) the 'k'
    // value of ionisation coeffs is shifted down  by one relative to the
    // recombination coeffs - therefore this
    // evaluation actually gives the balance

    iz_stage_distribution[k + 1] =
        iz_stage_distribution[k] * (k_iz_evaluated / k_rec_evaluated);

    sum_iz += iz_stage_distribution[k + 1];
  }

  // # Normalise such that the sum over all ionisation stages is '1' at all
  // points
  for (int k = 0; k <= Z; ++k) {
    iz_stage_distribution[k] = iz_stage_distribution[k] / sum_iz;
  }

  set<string> radiative_processes = {"line_power", "continuum_power"};
  if (impurity.get_has_charge_exchange()) {
    radiative_processes.insert("cx_power");
  }

  double total_power = 0;

  for (int k = 0; k < Z; ++k) {
    double k_power = 0;
    for (set<string>::iterator iter = radiative_processes.begin();
         iter != radiative_processes.end(); ++iter) {

      shared_ptr<RateCoefficient> rate_coefficient =
          impurity.get_rate_coefficient(*iter);
      double k_evaluated = rate_coefficient->call0D(k, Te, Ne);

      double scale;
      int target_charge_state;

      if (*iter == "line_power") {
        //# range of k is 0 to (Z-1)+ (needs bound electrons)
        target_charge_state = k; //#electron-bound target
        //# Prad = L * Ne * Nz^k+
        //#      = L * scale
        // N.b. Ne is function input
        double Nz_charge_state =
            Ni * iz_stage_distribution[target_charge_state];
        scale = Ne * Nz_charge_state;

      } else if (*iter == "continuum_power") {
        //# range of k is 1+ to Z+ (needs charged target)
        target_charge_state = k + 1; //#charged target
        //# Prad = L * Ne * Nz^(k+1)
        //#      = L * scale
        // N.b. Ne is function input
        double Nz_charge_state =
            Ni * iz_stage_distribution[target_charge_state];
        scale = Ne * Nz_charge_state;
      } else if (*iter == "cx_power") {
        //# range of k is 1+ to Z+ (needs charged target)
        target_charge_state = k + 1; //#charged target
        //# Prad = L * n_0 * Nz^(k+1)+
        //#      = L * scale
        // N.b. Nn is function input
        double Nz_charge_state =
            Ni * iz_stage_distribution[target_charge_state];
        scale = Nn * Nz_charge_state;
      } else {
        throw invalid_argument(
            "radiative_process not recognised (in computeRadiatedPower)");
      }
      double power = scale * k_evaluated;
      // N.b. These won't quite give the power from the kth charge state.
      // Instead they give the
      // power from the kth element on the rate coefficient, which may be kth or
      // (k+1)th charge state
      k_power += power;
    }
    // N.b. These won't quite give the power from the kth charge state. Instead
    // they give the
    // power from the kth element on the rate coefficient, which may be kth or
    // (k+1)th charge state
    // cout << "Total power from all procs. from k="<<k<<" is "<<k_power<<"
    // [W/m3]\n"<<endl;
    total_power += k_power;
  }

  return total_power;
}
