#include "../include/hydrogen_charge_exchange.hxx"

void HydrogenChargeExchange::calculate_rates(Options& atom1, Options& ion1,
                                             Options& atom2, Options& ion2,
                                             Field3D &R,
                                             Field3D &atom_mom, Field3D &ion_mom,
                                             Field3D &atom_energy, Field3D &ion_energy,
                                             Field3D &atom_rate, Field3D &ion_rate,
                                             BoutReal &rate_multiplier,
                                             bool &no_neutral_cx_mom_gain) {

  // Temperatures and masses of initial atom and ion
  const Field3D Tatom = GET_VALUE(Field3D, atom1["temperature"]);
  const BoutReal Aatom = GET_VALUE(BoutReal, atom1["AA"]);
  ASSERT1(GET_VALUE(BoutReal, ion2["AA"]) == Aatom); // Check that the mass is consistent

  const Field3D Tion = GET_VALUE(Field3D, ion1["temperature"]);
  const BoutReal Aion = GET_VALUE(BoutReal, ion1["AA"]);
  ASSERT1(GET_VALUE(BoutReal, atom2["AA"]) == Aion); // Check that the mass is consistent

  // Calculate effective temperature in eV
  Field3D Teff = (Tatom / Aatom + Tion / Aion) * Tnorm;
  for (auto& i : Teff.getRegion("RGN_NOBNDRY")) {
    if (Teff[i] < 0.01) {
      Teff[i] = 0.01;
    } else if (Teff[i] > 10000) {
      Teff[i] = 10000;
    }
  }
  const Field3D lnT = log(Teff);

  Field3D ln_sigmav = -18.5028;
  Field3D lnT_n = lnT; // (lnT)^n
  // b0 -1.850280000000E+01 b1 3.708409000000E-01 b2 7.949876000000E-03
  // b3 -6.143769000000E-04 b4 -4.698969000000E-04 b5 -4.096807000000E-04
  // b6 1.440382000000E-04 b7 -1.514243000000E-05 b8 5.122435000000E-07
  for (BoutReal b : {0.3708409, 7.949876e-3, -6.143769e-4, -4.698969e-4, -4.096807e-4,
                     1.440382e-4, -1.514243e-5, 5.122435e-7}) {
    ln_sigmav += b * lnT_n;
    lnT_n *= lnT;
  }

  // Get rate coefficient, convert cm^3/s to m^3/s then normalise
  // Optionally multiply by arbitrary multiplier
  const Field3D sigmav = exp(ln_sigmav) * (1e-6 * Nnorm / FreqNorm) * rate_multiplier;

  const Field3D Natom = floor(GET_VALUE(Field3D, atom1["density"]), 1e-5);
  const Field3D Nion = floor(GET_VALUE(Field3D, ion1["density"]), 1e-5);

  R = Natom * Nion * sigmav; // Rate coefficient. This is an output parameter.

  if ((&atom1 != &atom2) or (&ion1 != &ion2)) {
    // Transfer particles atom1 -> ion2, ion1 -> atom2
    subtract(atom1["density_source"], R);
    add(ion2["density_source"], R);
    subtract(ion1["density_source"], R);
    add(atom2["density_source"], R);
  } // Skip the case where the same isotope swaps places

  // Transfer momentum
  auto atom1_velocity = GET_NOBOUNDARY(Field3D, atom1["velocity"]);
  auto ion1_velocity = GET_NOBOUNDARY(Field3D, ion1["velocity"]);

  // Transfer fom atom1 to ion2
  atom_mom = R * Aatom * atom1_velocity;
  subtract(atom1["momentum_source"], atom_mom);
  if (no_neutral_cx_mom_gain == false) {
    add(ion2["momentum_source"], atom_mom);
  }

  // Transfer from ion1 to atom2
  ion_mom = R * Aion * ion1_velocity;
  subtract(ion1["momentum_source"], ion_mom);
  add(atom2["momentum_source"], ion_mom);

  // Frictional heating: Friction force between ions and atoms
  // converts kinetic energy to thermal energy
  //
  // This handles the general case that ion1 != ion2
  // and atom1 != atom2

  auto ion2_velocity = GET_NOBOUNDARY(Field3D, ion2["velocity"]);
  add(ion2["energy_source"], 0.5 * Aatom * R * SQ(ion2_velocity - atom1_velocity));

  auto atom2_velocity = GET_NOBOUNDARY(Field3D, atom2["velocity"]);
  add(atom2["energy_source"], 0.5 * Aion * R * SQ(atom2_velocity - ion1_velocity));

  // Transfer thermal energy
  atom_energy = (3. / 2) * R * Tatom;
  subtract(atom1["energy_source"], atom_energy);
  add(ion2["energy_source"], atom_energy);

  ion_energy = (3. / 2) * R * Tion;
  subtract(ion1["energy_source"], ion_energy);
  add(atom2["energy_source"], ion_energy);

  // Update collision frequency for the two colliding species in s^-1
  atom_rate = Nion * sigmav;
  ion_rate = Natom * sigmav;
  add(atom1["collision_frequency"], atom_rate);
  add(ion1["collision_frequency"], ion_rate);
}
