#include "../include/hydrogen_charge_exchange.hxx"

void HydrogenChargeExchange::calculate_rates(Options& atom1, Options& ion1,
                                             Options& atom2, Options& ion2) {

  // Temperatures and masses of initial atom and ion
  const Field3D Tatom = get<Field3D>(atom1["temperature"]);
  const BoutReal Aatom = get<BoutReal>(atom1["AA"]);
  ASSERT1(get<BoutReal>(ion2["AA"]) == Aatom); // Check that the mass is consistent

  const Field3D Tion = get<Field3D>(ion1["temperature"]);
  const BoutReal Aion = get<BoutReal>(ion1["AA"]);
  ASSERT1(get<BoutReal>(atom2["AA"]) == Aion); // Check that the mass is consistent

  // Calculate effective temperature in eV
  const Field3D Teff = (Tatom / Aatom + Tion / Aion) * Tnorm;
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
  const Field3D sigmav = exp(ln_sigmav) * (1e-6 * Nnorm / FreqNorm);

  const Field3D Natom = floor(get<Field3D>(atom1["density"]), 1e-5);
  const Field3D Nion = floor(get<Field3D>(ion1["density"]), 1e-5);

  const Field3D R = Natom * Nion * sigmav; // Rate coefficient

  if ((&atom1 != &atom2) or (&ion1 != &ion2)) {
    // Transfer particles atom1 -> ion2, ion1 -> atom2
    subtract(atom1["density_source"], R);
    add(ion2["density_source"], R);
    subtract(ion1["density_source"], R);
    add(atom2["density_source"], R);
  } // Skip the case where the same isotope swaps places

  // Transfer momentum
  const Field3D atom_mom = R * get<Field3D>(atom1["velocity"]);
  subtract(atom1["momentum_source"], atom_mom);
  add(ion2["momentum_source"], atom_mom);

  const Field3D ion_mom = R * get<Field3D>(ion1["velocity"]);
  subtract(ion1["momentum_source"], ion_mom);
  add(atom2["momentum_source"], ion_mom);

  // Transfer energy
  const Field3D atom_energy = (3. / 2) * R * Tatom;
  subtract(atom1["energy_source"], atom_energy);
  add(ion2["energy_source"], atom_energy);

  const Field3D ion_energy = (3. / 2) * R * Tion;
  subtract(ion1["energy_source"], ion_energy);
  add(atom2["energy_source"], ion_energy);
}
