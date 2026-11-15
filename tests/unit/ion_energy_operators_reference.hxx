#pragma once
// Independent per-cell/per-component reference formulas retained for regression.
#include "../../include/evolve_vorticity_ion_energy.hxx"
#include <bout/derivs.hxx>
#include <bout/region.hxx>
#include <algorithm>
#include <cmath>
#include <utility>
namespace ion_energy_reference {
struct CellFaces {
  BoutReal xp, xm, zp, zm;
};
inline CellFaces geometricPerpendicularFaceFluxes(const Field3D& a, const Field3D& f,
                                                  const Ind3D& i) {
  const auto* coord = f.getCoordinates();
  const auto& J = coord->J();
  const auto& g11 = coord->g11();
  const auto& g33 = coord->g33();
  const auto& g13 = coord->g13();
  const auto& dx = coord->dx();
  const auto& dz = coord->dz();

  const auto xflux = [&](const Ind3D& l) {
    const auto r = l.xp();
    return 0.5 * (a[l] + a[r]) * ((J[l] * g11[l] * dz[l]) + (J[r] * g11[r] * dz[r]))
           * (f[r] - f[l]) / (dx[l] + dx[r]);
  };
  const auto zflux = [&](const Ind3D& l) {
    const auto r = l.zp();
    return 0.5 * (a[l] + a[r]) * ((J[l] * g33[l] * dx[l]) + (J[r] * g33[r] * dx[r]))
           * (f[r] - f[l]) / (dz[l] + dz[r]);
  };
  const auto corner = [&](const Ind3D& l) {
    const auto x = l.xp();
    const auto z = l.zp();
    const auto xz = x.zp();
    const BoutReal coefficient =
        0.0625 * (a[l] + a[x] + a[z] + a[xz])
        * ((J[l] * g13[l]) + (J[x] * g13[x]) + (J[z] * g13[z]) + (J[xz] * g13[xz]));
    // Centre-to-centre lengths cancel against the corner area in the
    // integrated cross flux. These are averaged differences, not derivatives.
    return std::pair<BoutReal, BoutReal>{0.5 * coefficient * (f[x] - f[l] + f[xz] - f[z]),
                                         0.5 * coefficient
                                             * (f[z] - f[l] + f[xz] - f[x])};
  };
  const auto pp = corner(i);
  const auto mp = corner(i.xm());
  const auto pm = corner(i.zm());
  const auto mm = corner(i.xm().zm());
  return {xflux(i) + 0.5 * (pp.second + pm.second),
          xflux(i.xm()) + 0.5 * (mp.second + mm.second),
          zflux(i) + 0.5 * (pp.first + mp.first),
          zflux(i.zm()) + 0.5 * (pm.first + mm.first)};
}
inline CellFaces perpendicularFaceFluxes(const Field3D& a, const Field3D& f,
                                         const Ind3D& i, bool bndry_flux) {
  auto flux = geometricPerpendicularFaceFluxes(a, f, i);
  if (!bndry_flux && !f.getMesh()->periodicX) {
    if (f.getMesh()->firstX() && i.x() == f.getMesh()->xstart) {
      flux.xm = 0.0;
    }
    if (f.getMesh()->lastX() && i.x() == f.getMesh()->xend) {
      flux.xp = 0.0;
    }
  }
  return flux;
}

inline Field3D Div_a_Grad_perp_XZ(const Field3D& a, const Field3D& f, bool bndry_flux) {
  ASSERT0(a.getLocation() == f.getLocation());
  Field3D result{zeroFrom(f)};
  const auto* coord = f.getCoordinates();
  BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
    const auto flux = perpendicularFaceFluxes(a, f, i, bndry_flux);
    // A zero X derivative alone does not eliminate the g13 normal flux.
    result[i] = (flux.xp - flux.xm + flux.zp - flux.zm)
                / (coord->J()[i] * coord->dx()[i] * coord->dz()[i]);
  }
  return result;
}

inline Field3D kineticEnergyPerp(const Field3D& a, const Field3D& f) {
  ASSERT0(a.getLocation() == f.getLocation());
  Field3D result{zeroFrom(f)};
  const auto* coord = f.getCoordinates();
  BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
    const auto flux = perpendicularFaceFluxes(a, f, i, true);
    // Half of each face's 1/2 * delta(f) * flux belongs to this cell.
    result[i] = ((f[i.xp()] - f[i]) * flux.xp + (f[i] - f[i.xm()]) * flux.xm
                 + (f[i.zp()] - f[i]) * flux.zp + (f[i] - f[i.zm()]) * flux.zm)
                / (4. * coord->J()[i] * coord->dx()[i] * coord->dz()[i]);
  }
  return result;
}

} // namespace ion_energy_reference
