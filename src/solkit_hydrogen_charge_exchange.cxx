#include "../include/solkit_hydrogen_charge_exchange.hxx"

#include "../include/integrate.hxx"  // for cellAverage

void SOLKITHydrogenChargeExchange::calculate_rates(Options& atom, Options& ion) {
  const auto AA = get<BoutReal>(ion["AA"]);
  // Check that mass is consistent
  ASSERT1(get<BoutReal>(atom["AA"]) == AA);

  // Densities
  const Field3D Natom = floor(get<Field3D>(atom["density"]), 1e-5);
  const Field3D Nion = floor(get<Field3D>(ion["density"]), 1e-5);

  const Field3D Vatom = IS_SET(atom["velocity"]) ? GET_VALUE(Field3D, atom["velocity"]) : 0.0;
  const Field3D Vion = IS_SET(ion["velocity"]) ? GET_VALUE(Field3D, ion["velocity"]) : 0.0;

  const Field3D friction = cellAverage(
       [&](BoutReal natom, BoutReal nion, BoutReal vatom, BoutReal vion){
         // CONSTANT CROSS-SECTION 3E-19m2, COLD ION/NEUTRAL AND STATIC NEUTRAL ASSUMPTION
         auto R = natom * nion * 3e-19 * fabs(vion) *
           (Nnorm * rho_s0);
         return AA * (vion - vatom) * R;
       },
       Natom.getRegion("RGN_NOBNDRY"))(Natom, Nion, Vatom, Vion);

  // Transfer momentum
  subtract(ion["momentum_source"], friction);
  add(atom["momentum_source"], friction);
}
