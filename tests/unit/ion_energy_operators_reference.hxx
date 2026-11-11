#pragma once
// Independent per-cell/per-component reference formulas retained for regression.
#include "../../include/evolve_vorticity_ion_energy.hxx"
#include <bout/derivs.hxx>
#include <bout/region.hxx>
#include <algorithm>
#include <cmath>
#include <utility>
namespace ion_energy_reference {
using bout::ion_energy::XZFaceFluxes;
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
// Corner (i+1/2,k+1/2) stream function; two-point differences match XPPM.
inline BoutReal streamCorner(const Field3D& f, const Ind3D& i) {
  return 0.25 * (f[i] + f[i.xp()] + f[i.zp()] + f[i.xp().zp()]);
}

inline XZFaceFluxes advectionFaces(const Field3D& n, const Field3D& f, bool limited,
                                   bool close_radial_faces = true) {
  auto* m = f.getMesh();
  const auto* coord = f.getCoordinates();
  XZFaceFluxes flux{zeroFrom(f), zeroFrom(f)};
  const auto edge = [&](const Ind3D& i, bool x, bool right) {
    const BoutReal a = 2.0 * (n[x ? i.xp() : i.zp()] - n[i]);
    const BoutReal b = 0.5 * (n[x ? i.xp() : i.zp()] - n[x ? i.xm() : i.zm()]);
    const BoutReal c = 2.0 * (n[i] - n[x ? i.xm() : i.zm()]);
    const BoutReal slope =
        (a * b <= 0.0 || a * c <= 0.0)
            ? 0.0
            : std::copysign(std::min({std::abs(a), std::abs(b), std::abs(c)}), a);
    return n[i] + (right ? 0.5 : -0.5) * slope;
  };
  for (int x = m->xstart - 1; x <= m->xend; ++x) {
    for (int y = m->ystart; y <= m->yend; ++y) {
      for (int z = 0; z < m->LocalNz; ++z) {
        const Ind3D i =
            Ind3D((x * m->LocalNy + y) * m->LocalNz + z, m->LocalNy, m->LocalNz);
        const auto r = i.xp();
        if (!(close_radial_faces
              && ((x == m->xstart - 1 && m->firstX()) || (x == m->xend && m->lastX())))) {
          const BoutReal velocity = 0.5 * (coord->J()[i] + coord->J()[r])
                                    * (streamCorner(f, i) - streamCorner(f, i.zm()))
                                    / coord->dz()[i];
          const BoutReal nf =
              limited ? (velocity > 0 ? edge(i, true, true) : edge(r, true, false))
                      : 0.5 * (n[i] + n[r]);
          flux.x[i] = velocity * nf * coord->dy()[i] * coord->dz()[i];
        }
        if (x >= m->xstart) {
          const BoutReal velocity = coord->J()[i]
                                    * (streamCorner(f, i.xm()) - streamCorner(f, i))
                                    / coord->dx()[i];
          const BoutReal nf =
              limited ? (velocity > 0 ? edge(i, false, true) : edge(i.zp(), false, false))
                      : 0.5 * (n[i] + n[i.zp()]);
          flux.z[i] = velocity * nf * coord->dx()[i] * coord->dy()[i];
        }
      }
    }
  }
  return flux;
}

inline XZFaceFluxes solenoidalReynoldsFaceFluxes(const Field3D& c, const Field3D& phi,
                                                 const Field3D& Phi) {
  auto* m = phi.getMesh();
  const auto* coord = phi.getCoordinates();
  XZFaceFluxes flux{zeroFrom(phi), zeroFrom(phi)};
  Field3D stream = zeroFrom(phi);
  const auto average = [](const auto& f, const Ind3D& i) {
    return 0.25 * (f[i] + f[i.xp()] + f[i.zp()] + f[i.xp().zp()]);
  };
  const Field3D cJ = c * coord->J();
  for (int x = m->xstart - 1; x <= m->xend; ++x) {
    for (int y = m->ystart; y <= m->yend; ++y) {
      for (int z = 0; z < m->LocalNz; ++z) {
        const Ind3D i =
            Ind3D((x * m->LocalNy + y) * m->LocalNz + z, m->LocalNy, m->LocalNz);
        const auto gradient = [&](const Field3D& f) {
          return std::pair<BoutReal, BoutReal>{
              (f[i.xp()] - f[i] + f[i.xp().zp()] - f[i.zp()])
                  / (2 * average(coord->dx(), i)),
              (f[i.zp()] - f[i] + f[i.xp().zp()] - f[i.xp()])
                  / (2 * average(coord->dz(), i))};
        };
        const auto p = gradient(phi), h = gradient(Phi);
        stream[i] =
            -0.5 * average(cJ, i)
            * (average(coord->g11(), i) * p.first * h.first
               + average(coord->g13(), i) * (p.first * h.second + p.second * h.first)
               + average(coord->g33(), i) * p.second * h.second);
      }
    }
  }
  for (int x = m->xstart - 1; x <= m->xend; ++x) {
    for (int y = m->ystart; y <= m->yend; ++y) {
      for (int z = 0; z < m->LocalNz; ++z) {
        const Ind3D i =
            Ind3D((x * m->LocalNy + y) * m->LocalNz + z, m->LocalNy, m->LocalNz);
        flux.x[i] = (stream[i] - stream[i.zm()]) * coord->dy()[i];
        if (x >= m->xstart) {
          flux.z[i] -= (stream[i] - stream[i.xm()]) * coord->dy()[i];
        }
      }
    }
  }
  return flux;
}

inline XZFaceFluxes reynoldsFaces(const Field3D& c, const Field3D& vort,
                                  const Field3D& phi, const Field3D& psi,
                                  bool lambda_limited = true) {
  auto* m = phi.getMesh();
  const auto* coord = phi.getCoordinates();
  const Field3D Phi = phi + psi;
  const Field3D D{
      Coordinates::FieldMetric{coord->g11() * coord->g33() - SQ(coord->g13())}};
  Field3D lambda = 0.5 * c * D * coord->J()
                   * Div_a_Grad_perp_XZ(Field3D{1.0 / (coord->J() * D)}, phi, true);
  lambda.applyBoundary("neumann_o2");
  Field3D Q = bracket(phi, psi, BRACKET_ARAKAWA);
  Q.applyBoundary("neumann_o2");
  bout::ion_energy::applyVorticityWallBoundary(Q, true);
  if (!m->periodicX) {
    for (int x : {m->xstart, m->xend}) {
      const int step = x == m->xstart ? -1 : 1;
      const bool owns = step < 0 ? m->firstX() : m->lastX();
      if (!owns) {
        continue;
      }
      for (int y = m->ystart; y <= m->yend; ++y) {
        for (int z = 0; z < m->LocalNz; ++z) {
          for (int layer = 1; layer <= m->xstart; ++layer) {
            lambda(x + step * layer, y, z) =
                lambda(x, y, z) + layer * (lambda(x, y, z) - lambda(x - step, y, z));
          }
        }
      }
    }
  }
  m->communicate(lambda, Q);
  const auto omega =
      advectionFaces(0.5 * vort, phi, true, lambda_limited || !m->periodicX);
  const auto rotation =
      advectionFaces(lambda, Phi, lambda_limited, lambda_limited || !m->periodicX);
  const auto curl = solenoidalReynoldsFaceFluxes(c, phi, Phi);
  XZFaceFluxes result{omega.x + rotation.x + curl.x, omega.z + rotation.z + curl.z};
  const Field3D half_c = 0.5 * c;
  for (const auto& i : phi.getRegion("RGN_NOBNDRY")) {
    const auto q = perpendicularFaceFluxes(half_c, Q, i, false);
    result.x[i] += q.xp * coord->dy()[i];
    if (i.x() == m->xstart) {
      result.x[i.xm()] += q.xm * coord->dy()[i];
    }
    result.z[i] += q.zp * coord->dy()[i];
  }
  return result;
}
} // namespace ion_energy_reference
