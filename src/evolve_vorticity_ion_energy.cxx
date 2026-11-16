#include <bout/assert.hxx>
#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/field.hxx>
#include <bout/field3d.hxx>
#include <bout/fv_ops.hxx>
#include <bout/fv_ops_impl.hxx>
#include <bout/globals.hxx>
#include <bout/initialprofiles.hxx>
#include <bout/invert_pardiv.hxx>
#include <bout/options.hxx>
#include <bout/output.hxx>
#include <bout/output_bout_types.hxx>
#include <bout/region.hxx>
#include <bout/solver.hxx>
#include <bout/sys/range.hxx>
#include <bout/utils.hxx>
#include <bout/vecops.hxx>

#include "../include/component.hxx"
#include "../include/div_ops.hxx"
#include "../include/evolve_vorticity_ion_energy.hxx"
#include "../include/guarded_options.hxx"
#include "../include/hermes_build_config.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/permissions.hxx"

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

using bout::globals::mesh;

namespace {
/// Update every physical radial guard column, preserving interior and Y guards.
///
/// @param[in,out] field Field whose owned X-wall guards are overwritten at physical
///                     Y rows and every local Z index, using the field's mesh.
/// @param value Callback returning each guard value in field's units. Called as
///              value(interior_x, step, layer, y, z), where step points outwards
///              (-1 or +1) and layer starts at 1; interior_x is the wall-adjacent cell.
///
/// Returns no value. Periodic X is a no-op. No communication occurs; any further
/// side effects come from the callback.
template <typename Field, typename Value>
void fillRadialGuards(Field& field, Value value) {
  auto* localmesh = field.getMesh();
  if (localmesh->periodicX) {
    return;
  }
  const auto fill = [&](int interior, int step, int count) {
    for (int y = localmesh->ystart; y <= localmesh->yend; ++y) {
      for (int z = 0; z < localmesh->LocalNz; ++z) {
        for (int layer = 1; layer <= count; ++layer) {
          field(interior + step * layer, y, z) = value(interior, step, layer, y, z);
        }
      }
    }
  };
  if (localmesh->firstX()) {
    fill(localmesh->xstart, -1, localmesh->xstart);
  }
  if (localmesh->lastX()) {
    fill(localmesh->xend, 1, localmesh->LocalNx - localmesh->xend - 1);
  }
}

/// Average a local X-Y line over Z for wall-state and geometry calculations.
///
/// @param field Read-only potential, storage field or 2D/3D metric field.
/// @param x Local radial cell index.
/// @param y Local poloidal cell index.
/// @returns Arithmetic mean over all LocalNz cells, in the input field's units.
///          A 2D field is treated as constant in Z. Uses the global mesh and
///          assumes a complete local Z line; no communication or mutation occurs.
template <typename Field>
BoutReal meanZ(const Field& field, int x, int y) {
  BoutReal mean = 0.0;
  for (int z = 0; z < mesh->LocalNz; ++z) {
    mean += field(x, y, z);
  }
  return mean / mesh->LocalNz;
}

/// Average per-Y wall values over the connected flux surface with J dy weights.
///
/// @param values Read-only wall values indexed by local Y, with LocalNy entries.
///               Only physical Y rows contribute to the average.
/// @param interior Local index of the first interior radial cell at the wall.
/// @param guard Local index of its adjacent radial guard cell.
/// @returns Weighted surface mean in the input values' units, identical on all
///          ranks in mesh->getYcomm(interior). The face weight is the mean of
///          guard/interior J times interior dy, using Z-averaged geometry.
///
/// This MPI collective must be called by every rank in that Y communicator.
/// Callers restrict it to a first-X core surface and assume positive weights.
/// Uses the global mesh; values and geometry are preserved.
BoutReal fluxSurfaceAverage(const std::vector<BoutReal>& values, int interior,
                            int guard) {
  const auto* coord = mesh->getCoordinates();
  BoutReal local[2]{0.0, 0.0};
  for (int y = mesh->ystart; y <= mesh->yend; ++y) {
    const BoutReal weight =
        0.5 * (meanZ(coord->J(), guard, y) + meanZ(coord->J(), interior, y))
        * meanZ(coord->dy(), interior, y);
    local[0] += weight * values[static_cast<size_t>(y)];
    local[1] += weight;
  }
  BoutReal global[2];
  MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, mesh->getYcomm(interior));
  return global[0] / global[1];
}

/// Broadcast per-Y wall values into all physical radial guard columns and Z cells.
///
/// @param[in,out] field Destination storage, RHS or diagnostic field. Only radial
///                     guards at physical Y rows on this rank's X walls are written;
///                     interior cells and Y boundary cells retain their values.
/// @param inner Read-only inner wall values in field's units, indexed by local Y.
///              Must contain LocalNy entries when this rank owns the inner wall.
/// @param outer Read-only outer wall values, with the same convention as inner.
///
/// Uses the global mesh without communication. Periodic X leaves field unchanged.
/// Returns no value; the output is written directly into field.
void writeRadialValues(Field3D& field, const std::vector<BoutReal>& inner,
                       const std::vector<BoutReal>& outer) {
  if (mesh->periodicX) {
    return;
  }
  const auto write = [&](int begin, int end, const std::vector<BoutReal>& values) {
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      for (int x = begin; x < end; ++x) {
        for (int z = 0; z < mesh->LocalNz; ++z) {
          field(x, y, z) = values[static_cast<size_t>(y)];
        }
      }
    }
  };
  if (mesh->firstX()) {
    write(0, mesh->xstart, inner);
  }
  if (mesh->lastX()) {
    write(mesh->xend + 1, mesh->LocalNx, outer);
  }
}

/// Visit charged non-electron species through permission-checked option handles.
///
/// @param allspecies Species section; children named "e", without a charge, or
///                   with absolute charge below 1e-5 are skipped.
/// @param func Callback receiving each selected GuardedOptions child by value.
///             It may read or modify species state within its granted permissions.
///
/// Returns no value. Traversal reads options without writing fields or communicating;
/// state changes and other side effects are those of func.
template <typename Function>
void forAllIons(GuardedOptions allspecies, Function func) {
  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons
    }
    auto species = allspecies[kv.first];

    if (!IS_SET(species["charge"])) {
      continue;
    }
    if (fabs(get<BoutReal>(species["charge"])) < 1e-5) {
      // No charge
      continue;
    }
    func(species);
  }
}
/// Visit charged non-electron species in a read-only options tree.
///
/// @param allspecies Species section, with the same selection as the guarded overload:
///                   skip "e", missing charge, and absolute charge below 1e-5.
/// @param func Callback receiving each selected child as a const Options reference.
///
/// Returns no value. The species tree is preserved and no communication occurs;
/// func may have side effects on other state it captures.
template <typename Function>
void forAllIons(const Options& allspecies, Function func) {
  for (const auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue; // Skip electrons
    }
    const auto& species = allspecies[kv.first];

    if (!IS_SET(species["charge"])) {
      continue;
    }
    if (fabs(get<BoutReal>(species["charge"])) < 1e-5) {
      // No charge
      continue;
    }
    func(species);
  }
}

} // namespace

namespace bout::ion_energy {
/// Impose mean Neumann and zero AC midpoint conditions on a wall pressure field.
///
/// @param[in,out] f Pressure or pressure-potential field. Each owned radial guard
///                 layer is set to 2*DC(f)(interior,y) - f(interior,y,z), using the
///                 wall-adjacent interior cell for every layer.
///
/// Returns no value. Only physical Y rows at physical X walls are written;
/// interior cells and Y guards are preserved. Periodic X is a no-op and no
/// communication occurs.
/// @throws BoutException If nonperiodic X does not have complete local periodic Z lines.
void applyWallPressureBoundary(Field3D& f) {
  auto* localmesh = f.getMesh();
  if (localmesh->periodicX) {
    return;
  }
  if (localmesh->getNZPE() != 1 || localmesh->zstart != 0
      || localmesh->zend + 1 != localmesh->LocalNz) {
    throw BoutException("radial_boundary requires complete periodic Z lines");
  }
  // Calculate the mean once for each wall row, then copy the first guard outward.
  const Field2D mean = DC(f);
  fillRadialGuards(
      f, [&](int x, int, int, int y, int z) { return 2.0 * mean(x, y) - f(x, y, z); });
}

/// Fill vorticity guards for stress-free or Neumann radial walls.
///
/// @param[in,out] w Vorticity-like field whose owned radial guards are overwritten.
/// @param midpoint_zero True copies minus the wall-adjacent interior value, giving
///                      a zero face midpoint; false copies it with the same sign.
///
/// Returns no value. Every radial layer uses the same adjacent interior cell;
/// only physical Y rows are written. Interior cells and Y guards are preserved,
/// periodic X is a no-op, and no communication occurs.
void applyVorticityWallBoundary(Field3D& w, bool midpoint_zero) {
  fillRadialGuards(w, [&](int x, int, int, int y, int z) {
    return midpoint_zero ? -w(x, y, z) : w(x, y, z);
  });
}

/// Recover the inward wall-face derivative from a stored potential midpoint.
///
/// @param phi_midpoint Arithmetic guard/interior midpoint, usually Z-averaged.
/// @param phi0 Potential of the wall-adjacent interior cell in the same units.
/// @param dx0 Positive x-coordinate width of that interior cell.
/// @returns 2*(phi0 - phi_midpoint)/dx0, directed into the domain at either wall.
///          This is a coordinate derivative; the normal metric is not included.
///
/// No state is changed and no communication occurs.
BoutReal wallInwardDerivative(BoutReal phi_midpoint, BoutReal phi0, BoutReal dx0) {
  // The stored midpoint is the arithmetic mean of a mirrored guard/interior
  // pair. The finite-volume solve determines its wall-face gradient.
  return 2.0 * (phi0 - phi_midpoint) / dx0;
}

namespace {
/// Get or create a mesh-owned region for X faces or X-Z corners.
///
/// @param[in,out] m Mesh supplying index dimensions and owning the cached region.
/// @param name Cache key; an existing region is reused without checking its bounds.
/// @param low Inclusive lower X index, used only when creating the region.
/// @param high Inclusive upper X index, used only when creating the region.
/// @returns A mesh-owned reference covering the requested X range, physical Y
///          rows and all local Z cells, valid for the lifetime of m.
///
/// First use adds to the mesh region registry; no communication occurs. Mesh
/// ownership prevents stale pointer caches when a test mesh is destroyed.
const Region<Ind3D>& faceRegion(Mesh* m, const std::string& name, int low, int high) {
  if (!m->hasRegion3D(name)) {
    m->addRegion3D(name, Region<Ind3D>{low, high, m->ystart, m->yend, 0, m->LocalNz - 1,
                                       m->LocalNy, m->LocalNz});
  }
  return m->getRegion3D(name);
}
/// Get the local X faces adjoining physical cells, including both end faces.
///
/// @param[in,out] m Mesh owning the region cache.
/// @returns Mesh-owned indices xstart-1 through xend, physical Y and all local Z.
///          Index i denotes face (i+1/2,j,k) or corner (i+1/2,j,k+1/2).
///
/// First use creates the ion_energy_xfaces region; no communication occurs.
const Region<Ind3D>& xFaces(Mesh* m) {
  return faceRegion(m, "ion_energy_xfaces", m->xstart - 1, m->xend);
}
/// Get only the physical radial wall faces owned by this rank.
///
/// @param[in,out] m Mesh supplying wall ownership and owning the region cache.
/// @returns Indices at xstart-1 for an owned inner wall and xend for an owned outer
///          wall, at physical Y and all local Z. Periodic X or no owned walls gives
///          an empty region; periodic X uses a static, rather than mesh-owned, region.
///
/// First nonperiodic use creates the ion_energy_walls region in m. No communication
/// occurs, and inter-rank faces are excluded.
const Region<Ind3D>& wallFaces(Mesh* m) {
  if (m->periodicX) {
    static const Region<Ind3D> empty;
    return empty;
  }
  if (!m->hasRegion3D("ion_energy_walls")) {
    Region<Ind3D>::RegionIndices indices;
    for (int x : {m->xstart - 1, m->xend}) {
      if ((x == m->xstart - 1 && m->firstX()) || (x == m->xend && m->lastX())) {
        for (int y = m->ystart; y <= m->yend; ++y) {
          for (int z = 0; z < m->LocalNz; ++z) {
            indices.emplace_back((x * m->LocalNy + y) * m->LocalNz + z, m->LocalNy,
                                 m->LocalNz);
          }
        }
      }
    }
    m->addRegion3D("ion_energy_walls", Region<Ind3D>{indices});
  }
  return m->getRegion3D("ion_energy_walls");
}
/// Interpolate a cell-centred field to the upper X-Z corner with equal weights.
///
/// @param f Read-only 2D or 3D field with valid values at i, xp, zp and xp.zp.
/// @param i Lower-X, lower-Z cell of the four-cell corner stencil.
/// @returns Arithmetic corner mean in f's units; a 2D field is constant in Z.
///
/// Reads supplied guards as needed; no boundary filling, communication or mutation.
template <typename Field>
BoutReal cornerAverage(const Field& f, const Ind3D& i) {
  return 0.25 * (f[i] + f[i.xp()] + f[i.zp()] + f[i.xp().zp()]);
}
/// Construct the symmetric cross-metric contributions to a perpendicular flux.
///
/// @param a Read-only cell-centred coefficient, averaged to the corner.
/// @param f Read-only scalar whose transverse jumps supply the gradients.
/// @param coord Read-only geometry; J*g13 is averaged over the four corner cells.
/// @param i Lower-X, lower-Z cell of the corner stencil.
/// @returns Pair of corner contributions: first to the X flux from Z jumps,
///          second to the Z flux from X jumps. The transverse dy factor is omitted;
///          perpendicularX/Z average these corner values onto faces and apply dy.
///
/// Requires valid stencil values, including guards. No boundary filling,
/// communication or input mutation occurs.
std::pair<BoutReal, BoutReal> crossCorner(const Field3D& a, const Field3D& f,
                                          const Coordinates& coord, const Ind3D& i) {
  const auto x = i.xp(), z = i.zp(), xz = x.zp();
  const auto& J = coord.J();
  const auto& g13 = coord.g13();
  const BoutReal coefficient =
      0.25 * cornerAverage(a, i)
      * (J[i] * g13[i] + J[x] * g13[x] + J[z] * g13[z] + J[xz] * g13[xz]);
  return {0.5 * coefficient * (f[z] - f[i] + f[xz] - f[x]),
          0.5 * coefficient * (f[x] - f[i] + f[xz] - f[z])};
}
/// Area-integrated diagonal X flux, including dy and the physical wall jump.
/// Reads only the two adjoining cells; no boundary filling or communication.
BoutReal perpendicularXDiagonal(const Field3D& a, const Field3D& f,
                                const Coordinates& coord, const Ind3D& i) {
  const auto r = i.xp();
  return 0.5 * (a[i] + a[r])
         * (coord.J()[i] * coord.g11()[i] * coord.dz()[i]
            + coord.J()[r] * coord.g11()[r] * coord.dz()[r])
         * (f[r] - f[i]) / (coord.dx()[i] + coord.dx()[r]) * coord.dy()[i];
}
/// Area-integrated diagonal Z flux, including dy.
/// Reads only the two adjoining cells; no boundary filling or communication.
BoutReal perpendicularZDiagonal(const Field3D& a, const Field3D& f,
                                const Coordinates& coord, const Ind3D& i) {
  const auto r = i.zp();
  return 0.5 * (a[i] + a[r])
         * (coord.J()[i] * coord.g33()[i] * coord.dx()[i]
            + coord.J()[r] * coord.g33()[r] * coord.dx()[r])
         * (f[r] - f[i]) / (coord.dz()[i] + coord.dz()[r]) * coord.dy()[i];
}
/// Combine diagonal and corner cross terms on one increasing-X face.
///
/// @param a Read-only cell-centred coefficient of a*grad_perp(f).
/// @param f Read-only scalar with a valid jump from i to i.xp().
/// @param cross Read-only crossCorner contributions stored at corner indices;
///              cross.x at i and i.zm() must be defined, without a dy factor.
/// @param coord Read-only J, g11 and cell widths at the face stencil.
/// @param i Lower-X cell adjoining face (i+1/2,j,k).
/// @returns Area-integrated flux towards increasing X, including dy and g13.
///
/// No wall closure, boundary filling, communication or input mutation occurs.
BoutReal perpendicularX(const Field3D& a, const Field3D& f, const XZFaceFluxes& cross,
                        const Coordinates& coord, const Ind3D& i) {
  return perpendicularXDiagonal(a, f, coord, i)
         + 0.5 * (cross.x[i] + cross.x[i.zm()]) * coord.dy()[i];
}
/// Combine diagonal and corner cross terms on one increasing-Z face.
///
/// @param a Read-only cell-centred coefficient of a*grad_perp(f).
/// @param f Read-only scalar with a valid jump from i to i.zp().
/// @param cross Read-only crossCorner contributions stored at corner indices;
///              cross.z at i and i.xm() must be defined, without a dy factor.
/// @param coord Read-only J, g33 and cell widths at the face stencil.
/// @param i Lower-Z cell adjoining face (i,j,k+1/2).
/// @returns Area-integrated flux towards increasing Z, including dy and g13.
///
/// No wall closure, boundary filling, communication or input mutation occurs.
BoutReal perpendicularZ(const Field3D& a, const Field3D& f, const XZFaceFluxes& cross,
                        const Coordinates& coord, const Ind3D& i) {
  return perpendicularZDiagonal(a, f, coord, i)
         + 0.5 * (cross.z[i] + cross.z[i.xm()]) * coord.dy()[i];
}
/// Reconstruct the donor-side scalar at a face with the MC slope limiter.
///
/// @tparam direction Stencil direction, X or Z.
/// @param n Read-only cell-centred scalar with a valid three-point donor stencil.
/// @param i Lower-coordinate cell adjoining the requested face.
/// @param velocity Signed face velocity or integrated velocity; only its sign is
///                 used. Positive selects i's right edge; otherwise the upper
///                 neighbour's left edge is used, including when velocity is zero.
/// @returns Reconstructed scalar in n's units, before multiplication by velocity.
///
/// X-face reconstruction can read two radial guard layers. No boundary filling,
/// communication or input mutation occurs.
template <DIRECTION direction>
BoutReal upwind(const Field3D& n, const Ind3D& i, BoutReal velocity) {
  const auto donor = velocity > 0 ? i : i.plus<1, direction>();
  FV::Stencil1D stencil;
  stencil.m = n[donor.template minus<1, direction>()];
  stencil.c = n[donor];
  stencil.p = n[donor.template plus<1, direction>()];
  FV::MC{}(stencil);
  return velocity > 0 ? stencil.R : stencil.L;
}
} // namespace

/// Suppress a flux through physical X walls while retaining all other faces.
///
/// @param[in,out] flux Face fields in the XZFaceFluxes index convention. Only flux.x
///                    at owned physical wall rows and physical Y is set to zero;
///                    flux.z, interior faces and inter-rank faces are preserved.
///
/// Returns no value. Allocates or detaches flux.x storage even for periodic X,
/// whose values are unchanged. May create the mesh's wall-region cache; no
/// communication occurs.
void closeRadialFaces(XZFaceFluxes& flux) {
  flux.x.allocate();
  BOUT_FOR(i, wallFaces(flux.x.getMesh())) { flux.x[i] = 0.0; }
}
/// Construct the shared finite-volume faces for a*grad_perp(f), including g13.
///
/// @param a Read-only cell-centred coefficient at the same location as f.
/// @param f Read-only scalar with valid adjacent guard values and geometry.
/// @returns Area-integrated fluxes, including dy, towards increasing X and Z.
///          x(i,j,k) is face (i+1/2,j,k), defined from xstart-1 through xend;
///          z(i,j,k) is face (i,j,k+1/2), defined on physical cells. Both use
///          physical Y and all local Z; unused guard entries are zero.
///
/// Physical wall gradients are retained. No boundary filling, communication or
/// input mutation occurs; the mesh's face-region cache may be initialised.
XZFaceFluxes perpendicularFaces(const Field3D& a, const Field3D& f) {
  ASSERT0(a.getLocation() == f.getLocation());
  const auto* coord = f.getCoordinates();
  XZFaceFluxes cross{zeroFrom(f), zeroFrom(f)}, flux{zeroFrom(f), zeroFrom(f)};
  BOUT_FOR(i, xFaces(f.getMesh())) {
    const auto value = crossCorner(a, f, *coord, i);
    cross.x[i] = value.first;
    cross.z[i] = value.second;
  }
  BOUT_FOR(i, xFaces(f.getMesh())) { flux.x[i] = perpendicularX(a, f, cross, *coord, i); }
  BOUT_FOR(i, f.getRegion("RGN_NOBNDRY")) {
    flux.z[i] = perpendicularZ(a, f, cross, *coord, i);
  }
  return flux;
}
/// Evaluate div(a*grad_perp(f)) with the shared conservative X-Z face stencil.
///
/// @param a Read-only cell-centred coefficient at the same location as f.
/// @param f Read-only scalar with valid guards for perpendicularFaces.
/// @param bndry_flux True retains geometric physical X-wall fluxes; false sets
///                   only those faces to zero, preserving periodic and inter-rank faces.
/// @returns Divergence per cell volume on physical cells, with unused guards zero.
///
/// No input mutation, boundary filling or communication occurs. Temporary faces
/// may be closed, and mesh face/wall-region caches may be initialised.
Field3D Div_a_Grad_perp_XZ(const Field3D& a, const Field3D& f, bool bndry_flux) {
  auto flux = perpendicularFaces(a, f);
  if (!bndry_flux) {
    closeRadialFaces(flux);
  }
  return divergence(flux);
}
/// Calculate perpendicular kinetic energy consistent with the forward face operator.
///
/// @param a Read-only kinetic-energy coefficient, at the same location as f.
/// @param f Read-only total drift potential, with valid adjacent guard values.
/// @returns 0.5*faceWork(perpendicularFaces(a,f),f) on physical cells, with unused
///          guards zero. Each adjacent face contributes one quarter of its
///          potential jump times integrated flux, divided by cell volume.
///
/// Uses geometric faces, including wall gradients, without wall closure. No input
/// mutation, boundary filling or communication occurs; mesh face-region caches
/// may be initialised.
Field3D kineticEnergyPerp(const Field3D& a, const Field3D& f) {
  // faceWork assigns half the jump*flux to each cell; energy supplies another half.
  return 0.5 * faceWork(perpendicularFaces(a, f), f);
}
/// Take a conservative cell-volume divergence, optionally weighting each face.
///
/// @param flux Read-only area-integrated, increasing-coordinate face fluxes. Each
///             physical cell's two X faces and periodic Z faces must be defined.
/// @param weight Optional read-only cell-centred scalar, averaged arithmetically
///               across each face. Inter-rank guards must already be communicated;
///               physical radial guards are replaced by Neumann values in a copy.
/// @returns Net outward weighted flux divided by J*dx*dy*dz on physical cells,
///          with unused guards zero. nullptr uses the supplied flux without weighting.
///
/// Inputs are preserved. Only the temporary weight copy's guards are changed;
/// no communication or flux closure occurs.
Field3D divergence(const XZFaceFluxes& flux, const Field3D* weight) {
  Field3D result = zeroFrom(flux.x);
  const auto* coord = result.getCoordinates();
  if (weight == nullptr) {
    BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
      result[i] = (flux.x[i] - flux.x[i.xm()] + flux.z[i] - flux.z[i.zm()])
                  / (coord->J()[i] * coord->dx()[i] * coord->dy()[i] * coord->dz()[i]);
    }
  } else {
    Field3D w = copy(*weight);
    fillRadialGuards(w, [&](int x, int, int, int y, int z) { return w(x, y, z); });
    BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
      result[i] =
          0.5
          * ((w[i] + w[i.xp()]) * flux.x[i] - (w[i] + w[i.xm()]) * flux.x[i.xm()]
             + (w[i] + w[i.zp()]) * flux.z[i] - (w[i] + w[i.zm()]) * flux.z[i.zm()])
          / (coord->J()[i] * coord->dx()[i] * coord->dy()[i] * coord->dz()[i]);
    }
  }
  return result;
}
/// Allocate half of each face's flux-times-potential-jump work to each adjacent cell.
///
/// @param flux Read-only area-integrated fluxes towards increasing X and Z, with
///             both X faces and periodic Z faces of every physical cell defined.
/// @param f Read-only potential with valid neighbours, including wall guards.
/// @returns Half the sum of the four integrated face fluxes times their respective
///          increasing-coordinate jumps in f, divided by J*dx*dy*dz. Defined on
///          physical cells, with unused guards zero.
///
/// No input mutation, boundary filling, communication or flux closure occurs.
Field3D faceWork(const XZFaceFluxes& flux, const Field3D& f) {
  Field3D result = zeroFrom(f);
  const auto* coord = f.getCoordinates();
  BOUT_FOR(i, result.getRegion("RGN_NOBNDRY")) {
    result[i] = 0.5
                * ((f[i.xp()] - f[i]) * flux.x[i] + (f[i] - f[i.xm()]) * flux.x[i.xm()]
                   + (f[i.zp()] - f[i]) * flux.z[i] + (f[i] - f[i.zm()]) * flux.z[i.zm()])
                / (coord->J()[i] * coord->dx()[i] * coord->dy()[i] * coord->dz()[i]);
  }
  return result;
}
namespace {
/// Static coefficients of J*R after replacing pure tangential derivatives with
/// diagonal FV divergences. The order is Phi_x*phi_x, Phi_x*phi_z,
/// Phi_z*phi_x, Phi_z*phi_z. These expressions follow by symbolic collection
/// of G_impl-omega*v_E, including the diagonal FV metric-gradient correction.
template <typename Field>
std::array<Field, 8> reynoldsMetricCoefficients(const Field& c) {
  auto* m = c.getMesh();
  ASSERT0(m->xstart >= 2);
  const auto* coord = c.getCoordinates();
  std::array<Field, 8> coeff;
  for (auto& field : coeff) {
    field = zeroFrom(c);
  }
  const Field cJ = c * coord->J();
  const Coordinates::FieldMetric det = coord->g11() * coord->g33() - SQ(coord->g13());
  BOUT_FOR(
      i, faceRegion(m, "ion_energy_reynolds_coefficients", m->xstart - 1, m->xend + 1)) {
    if constexpr (std::is_same_v<Field, Field2D>) {
      if (i.z() != 0) {
        continue;
      }
    }
    const auto derivativeX = [&](const auto& f) {
      return (f[i.xp()] - f[i.xm()]) / (2 * coord->dx()[i]);
    };
    const auto derivativeZ = [&](const auto& f) {
      return (f[i.zp()] - f[i.zm()]) / (2 * coord->dz()[i]);
    };
    const BoutReal A = cJ[i], Ax = derivativeX(cJ), Az = derivativeZ(cJ);
    const BoutReal D = det[i], Dx = derivativeX(det), Dz = derivativeZ(det);
    const BoutReal g11 = coord->g11()[i], g13 = coord->g13()[i], g33 = coord->g33()[i];
    const BoutReal g11x = derivativeX(coord->g11()), g11z = derivativeZ(coord->g11());
    const BoutReal g13x = derivativeX(coord->g13()), g13z = derivativeZ(coord->g13());
    const BoutReal g33x = derivativeX(coord->g33()), g33z = derivativeZ(coord->g33());
    coeff[0][i] = -1.0 / 2.0 * A * g11z - 1.0 / 2.0 * Az * g11;
    coeff[1][i] = -1.0 / 2.0 * A * g11x - A * g13z - 1.0 / 2.0 * Ax * g11 - Az * g13;
    coeff[2][i] = (1.0 / 2.0)
                  * (A * D * g11x * g33 + 2 * A * D * g13 * g33z - A * Dx * g11 * g33
                     - A * Dz * g13 * g33 + Az * D * g13 * g33)
                  / (D * g33);
    coeff[3][i] =
        (1.0 / 2.0) * (A * D * g33z - A * Dx * g13 - A * Dz * g33 - Ax * D * g13) / D;
    coeff[4][i] =
        -1.0 / 2.0 * (A * D * g11x - A * Dx * g11 - A * Dz * g13 - Az * D * g13) / D;
    coeff[5][i] = -1.0 / 2.0
                  * (A * D * g11 * g33z + 2 * A * D * g11x * g13 - A * Dx * g11 * g13
                     - A * Dz * g11 * g33 + Ax * D * g11 * g13)
                  / (D * g11);
    coeff[6][i] = A * g13x + (1.0 / 2.0) * A * g33z + Ax * g13 + (1.0 / 2.0) * Az * g33;
    coeff[7][i] = (1.0 / 2.0) * A * g33x + (1.0 / 2.0) * Ax * g33;
  }
  return coeff;
}

/// Evaluate the shared current using precomputed, immutable metric coefficients.
template <typename Field>
XZFaceFluxes reynoldsFaces(const Field3D& c, const Field3D& vort, const Field3D& phi,
                           const Field3D& psi, const std::array<Field, 8>& coeff) {
  auto* m = phi.getMesh();
  ASSERT0(m->xstart >= 2);
  const auto* coord = phi.getCoordinates();
  const Field3D Phi = phi + psi;
  Field3D phi_corner = zeroFrom(phi), Phi_corner = zeroFrom(phi);
  XZFaceFluxes gradient{zeroFrom(phi), zeroFrom(phi)}, flux{zeroFrom(phi), zeroFrom(phi)};
  BOUT_FOR(i, xFaces(m)) {
    phi_corner[i] = cornerAverage(phi, i);
    Phi_corner[i] = cornerAverage(Phi, i);
    gradient.x[i] = (Phi[i.xp()] - Phi[i] + Phi[i.xp().zp()] - Phi[i.zp()])
                    / (2 * cornerAverage(coord->dx(), i));
    gradient.z[i] = (Phi[i.zp()] - Phi[i] + Phi[i.xp().zp()] - Phi[i.xp()])
                    / (2 * cornerAverage(coord->dz(), i));
  }
  const auto volume = [&](const Ind3D& i) {
    return coord->J()[i] * coord->dx()[i] * coord->dy()[i] * coord->dz()[i];
  };
  // Evaluate each diagonal divergence directly from original potential samples.
  // In particular Lx uses the physical wall-face jumps; no derived guards exist.
  const auto Lx = [&](const Ind3D& i) {
    return (perpendicularXDiagonal(c, Phi, *coord, i)
            - perpendicularXDiagonal(c, Phi, *coord, i.xm()))
           / volume(i);
  };
  const auto Lz = [&](const Ind3D& i) {
    return (perpendicularZDiagonal(c, Phi, *coord, i)
            - perpendicularZDiagonal(c, Phi, *coord, i.zm()))
           / volume(i);
  };
  const auto bilinear = [&](const auto& coefficients, size_t offset, const Ind3D& i,
                            const Ind3D& r, BoutReal hx, BoutReal hz, BoutReal px,
                            BoutReal pz) {
    return 0.5
           * ((coefficients[offset + 0][i] + coefficients[offset + 0][r]) * hx * px
              + (coefficients[offset + 1][i] + coefficients[offset + 1][r]) * hx * pz
              + (coefficients[offset + 2][i] + coefficients[offset + 2][r]) * hz * px
              + (coefficients[offset + 3][i] + coefficients[offset + 3][r]) * hz * pz);
  };
  BOUT_FOR(i, xFaces(m)) {
    const auto r = i.xp();
    const BoutReal J = 0.5 * (coord->J()[i] + coord->J()[r]);
    const BoutReal dz = 0.5 * (coord->dz()[i] + coord->dz()[r]);
    const BoutReal dx = 0.5 * (coord->dx()[i] + coord->dx()[r]);
    const BoutReal velocity = J * (phi_corner[i] - phi_corner[i.zm()]);
    const BoutReal px = (phi[r] - phi[i]) / dx, pz = velocity / (J * dz);
    const BoutReal hx = (Phi[r] - Phi[i]) / dx;
    const BoutReal hz = (Phi_corner[i] - Phi_corner[i.zm()]) / dz;
    const BoutReal hxz = (gradient.x[i] - gradient.x[i.zm()]) / dz;
    const BoutReal cf = 0.5 * (c[i] + c[r]);
    const BoutReal g11 = 0.5 * (coord->g11()[i] + coord->g11()[r]);
    const BoutReal g13 = 0.5 * (coord->g13()[i] + coord->g13()[r]);
    const BoutReal g33 = 0.5 * (coord->g33()[i] + coord->g33()[r]);
    const BoutReal tangential = 0.5 * (Lz(i) + Lz(r));
    const BoutReal R = -J * cf * (g11 * px + g13 * pz) * hxz
                       - J * (g13 / g33) * px * tangential
                       + bilinear(coeff, 0, i, r, hx, hz, px, pz);
    // The same velocity multiplies omega and the tangential FV diagonal part.
    // Central omega would cancel that part exactly, including at shared faces.
    flux.x[i] =
        (velocity * (upwind<DIRECTION::X>(vort, i, velocity) - tangential) + dz * R)
        * coord->dy()[i];
  }
  BOUT_FOR(i, phi.getRegion("RGN_NOBNDRY")) {
    const auto r = i.zp();
    // Retain the existing integrated Z velocity, including its cell-centred J.
    const BoutReal J = coord->J()[i], dx = coord->dx()[i];
    const BoutReal dz = 0.5 * (coord->dz()[i] + coord->dz()[r]);
    const BoutReal velocity = J * (phi_corner[i.xm()] - phi_corner[i]);
    const BoutReal px = -velocity / (J * dx), pz = (phi[r] - phi[i]) / dz;
    const BoutReal hx = (Phi_corner[i] - Phi_corner[i.xm()]) / dx;
    const BoutReal hz = (Phi[r] - Phi[i]) / dz;
    const BoutReal hxz = (gradient.z[i] - gradient.z[i.xm()]) / dx;
    const BoutReal cf = 0.5 * (c[i] + c[r]);
    const BoutReal g11 = 0.5 * (coord->g11()[i] + coord->g11()[r]);
    const BoutReal g13 = 0.5 * (coord->g13()[i] + coord->g13()[r]);
    const BoutReal g33 = 0.5 * (coord->g33()[i] + coord->g33()[r]);
    const BoutReal tangential = 0.5 * (Lx(i) + Lx(r));
    const BoutReal R = J * cf * (g13 * px + g33 * pz) * hxz
                       + J * (g13 / g11) * pz * tangential
                       + bilinear(coeff, 4, i, r, hx, hz, px, pz);
    flux.z[i] =
        (velocity * (upwind<DIRECTION::Z>(vort, i, velocity) - tangential) + dx * R)
        * coord->dy()[i];
  }
  return flux;
}
} // namespace

/// Construct full-vorticity MC advection plus the compact central remainder R.
///
/// @param c Static Boussinesq coefficient, with valid geometry and guards.
/// @param vort Full physical vorticity, with two valid radial MC guard layers.
/// @param phi Electric potential, with communicated guards and reflected walls.
/// @param psi Pressure potential, with communicated guards and reflected AC walls.
/// @returns Area-integrated increasing-coordinate faces including dy: X indices
///          xstart-1..xend, Z on physical cells; unused entries are zero.
///
/// R reads only one radial guard layer of the input potentials. Tangential
/// diagonal FV divergences are evaluated from those samples, never extrapolated
/// or communicated. The MC stencil of vort still reads two layers. Reflected
/// wall midpoints make the advective X-wall velocity zero; R retains the same
/// two-point normal wall gradients as perpendicularFaces and their Z differences.
/// No wall closure, input mutation, boundary filling or communication occurs.
/// This standalone entry computes geometry coefficients locally; the component
/// precomputes them once in its constructor and uses the same evaluator.
/// Coefficient gradients use two supplied geometry/coefficient guard layers;
/// no derived field guard values are filled or read by the current evaluator.
XZFaceFluxes reynoldsFaceFluxes(const Field3D& c, const Field3D& vort, const Field3D& phi,
                                const Field3D& psi) {
  return reynoldsFaces(c, vort, phi, psi, reynoldsMetricCoefficients(c));
}
} // namespace bout::ion_energy

using bout::ion_energy::Div_a_Grad_perp_XZ;

EvolveVorticityIonEnergy::EvolveVorticityIonEnergy(std::string name, Options& alloptions,
                                                   Solver* solver)
    : NamedComponent(
          name,
          {readOnly("species:{all_species}:charge"),
           readIfSet("species:{charged}:velocity", Regions::Interior),
           readIfSet("species:{charged}:momentum", Regions::Interior),

           // Electron inputs
           readOnly("species:e:pressure", Regions::Interior), readIfSet("species:e:AA"),
           // Electron outputs
           readWrite("species:e:energy_source"),
           // Ion inputs
           readOnly("species:{ions}:density", Regions::Interior),
           readOnly("species:{ions}:AA"),
           // Ion outputs
           readWrite("species:{ions}:pressure"), readWrite("species:{ions}:temperature"),
           // Fields
           readWrite("fields:vorticity"), readWrite("fields:phi"),
           readWrite("fields:ion_boussinesq_coefficient"),
           readWrite("fields:ion_pressure_potential"),
           readWrite("fields:ion_mass_density"), readWrite("fields:ion_reynolds_x"),
           readWrite("fields:ion_reynolds_z"), readWrite("fields:DivJReynolds"),
           readWrite("fields:DivJdia"), readWrite("fields:ion_wall_rate")}) {

  auto& options = alloptions[name];
  // Conservative drift is an explicit companion, not an implicit second RHS.
  int drift_count = 0;
  std::vector<std::string> enabled;
  const Options& scheduler =
      alloptions.isSet("components") ? alloptions : alloptions["hermes"];
  if (scheduler.isSet("components")) {
    for (const auto& item : strsplit(scheduler["components"].as<std::string>(), ',')) {
      enabled.push_back(trim(item, " \t\r\n()"));
    }
  }
  for (const auto& [component_name, section] : alloptions.getChildren()) {
    if (!enabled.empty()
        && std::find(enabled.begin(), enabled.end(), component_name) == enabled.end()) {
      continue;
    }
    const std::string component_type =
        section.isSet("type") ? section["type"].as<std::string>() : component_name;
    if (component_type == "polarisation_drift") {
      if (section["implementation"].withDefault(std::string{"legacy"})
          != "conservative") {
        throw BoutException("evolve_vorticity_ion_energy requires "
                            "polarisation_drift:implementation=conservative");
      }
      ++drift_count;
    }
  }
  if (drift_count != 1) {
    throw BoutException("evolve_vorticity_ion_energy requires exactly one enabled "
                        "conservative polarisation_drift component");
  }
  for (const auto* key : {"polarisation_ac_boundary", "polarisation_corrections"}) {
    if (options.isSet(key)) {
      throw BoutException("{}:{} moved to the conservative polarisation_drift section",
                          name, key);
    }
  }
  if (options.isSection("polarisation_laplacian")
      && !options["polarisation_laplacian"].getChildren().empty()) {
    throw BoutException("{}:polarisation_laplacian moved to polarisation_drift:laplacian",
                        name);
  }

  const Options& units = alloptions["units"];
  const BoutReal Tnorm = units["eV"];
  const BoutReal Omega_ci = 1. / units["seconds"].as<BoutReal>();
  const BoutReal Bnorm = units["Tesla"];
  const BoutReal Lnorm = units["meters"];

  using bout::ion_energy::PolarisationACBoundary;
  using bout::ion_energy::RadialBoundary;
  radial_model.walls = !mesh->periodicX;
  radial_model.mode = options["radial_boundary"]
                          .doc("Radial wall model: fixed, free_slip or no_slip")
                          .withDefault(RadialBoundary::fixed);
  radial_model.evolve_wall_potential = radial_model.mode == RadialBoundary::no_slip;
  radial_model.vorticity_wall_zero = radial_model.mode == RadialBoundary::free_slip;
  radial_model.viscous_wall_flux = radial_model.vorticity_wall_zero;
  radial_model.core_flux_surface_average =
      radial_model.walls && mesh->firstX() && mesh->periodicY(mesh->xstart);

  // These flags are owned by the wall model, including when X is periodic.
  // Reject before creating a solver, which would otherwise consume the setting.
  for (const auto* section : {"laplacian"}) {
    for (const auto* key : {"inner_boundary_flags", "outer_boundary_flags"}) {
      if (options[section][key].isSet()) {
        throw BoutException("{}:{} is controlled by radial_boundary", section, key);
      }
    }
  }
  auto& vort_options = Options::root()["Vort"];
  for (const auto* key : {"bndry_xin", "bndry_xout"}) {
    if (vort_options[key].isSet()) {
      throw BoutException("Vort:{} is controlled by radial_boundary", key);
    }
  }
  auto& evolve_boundary = vort_options["evolve_bndry"];
  if (evolve_boundary.isSet() && !evolve_boundary.as<bool>()) {
    throw BoutException("radial_boundary requires Vort:evolve_bndry=true");
  }
  if (radial_model.evolve_wall_potential) {
    const BoutReal kappa_SI =
        options["phi_boundary_velocity"]
            .doc("Required finite positive velocity relaxing the mean wall-normal "
                 "E-field to zero for radial_boundary=no_slip [m/s].")
            .as<BoutReal>();
    if (!std::isfinite(kappa_SI) || kappa_SI <= 0.0) {
      throw BoutException("radial_boundary=no_slip requires finite positive "
                          "phi_boundary_velocity");
    }
    radial_model.phi_boundary_velocity = kappa_SI / (Lnorm * Omega_ci);
  }
  if (radial_model.walls) {
    if (mesh->xstart < 2) {
      throw BoutException("radial_boundary requires two radial guard cells");
    }
    evolve_boundary.force(true, "evolve_vorticity_ion_energy");
    wall_rate_inner.resize(static_cast<size_t>(mesh->LocalNy), 0.0);
    wall_rate_outer.resize(static_cast<size_t>(mesh->LocalNy), 0.0);
    // Own a no-op boundary configuration before moving the original Vort
    // boundary operators into the physical-vorticity field. Swap only the
    // FieldData metadata: solver pointers to Vort and ddt(Vort) stay valid.
    const std::string storage_name = name + "_boundary_storage";
    Options::root()[storage_name]["bndry_all"].force(std::string{"none"});
    vorticity.setBoundary(storage_name);
  }
  output_info.write(
      "{}: radial_boundary={}, walls={}, "
      "evolve_wall_potential={}, vorticity_wall_zero={}, viscous_wall_flux={}, "
      "core_flux_surface_average={}, "
      "phi_boundary_velocity(normalised)={}, phi_flags={}\n",
      name, toString(radial_model.mode), radial_model.walls,
      radial_model.evolve_wall_potential, radial_model.vorticity_wall_zero,
      radial_model.viscous_wall_flux, radial_model.core_flux_surface_average,
      radial_model.phi_boundary_velocity, radial_model.phi_flags);

  // Evolve the vorticity and ion energy in time
  solver->add(Vort, "Vort");
  solver->add(Eion, "Eion");

  if (radial_model.walls) {
    swap(static_cast<FieldData&>(Vort), static_cast<FieldData&>(vorticity));
    ddt(Vort).copyBoundary(Vort);
    initialiseBoundaryState();
  }

  const auto* coord = mesh->getCoordinates();

  average_atomic_mass = options["average_atomic_mass"]
                            .doc("Fixed reference atomic mass for the Boussinesq "
                                 "polarisation coefficient")
                            .withDefault<BoutReal>(2.0); // Deuterium

  poloidal_flows =
      options["poloidal_flows"].doc("Include poloidal ExB flow").withDefault<bool>(false);
  phi_dissipation = options["phi_dissipation"]
                        .doc("Numerically damp parallel potential variations using "
                             "the completed collective sound_speed")
                        .withDefault<bool>(false);

  density_floor = options["density_floor"].doc("Minimum density floor").withDefault(1e-7);

  // Match vorticity's physical units, scalar coefficient and guards.
  include_viscosity = options.isSet("viscosity");
  viscosity =
      options["viscosity"].doc("Kinematic viscosity [m^2/s]").withDefault<BoutReal>(0.0)
      / (Lnorm * Lnorm * Omega_ci);
  viscosity.applyBoundary("dirichlet");
  // Keep the physical coefficient at radial faces. Dirichlet coefficient guards
  // would cancel nu at the wall and suppress even the free-slip viscous flux.
  if (radial_model.walls) {
    fillRadialGuards(viscosity,
                     [&](int x, int, int, int y, int z) { return viscosity(x, y, z); });
  }
  mesh->communicate(viscosity);
  if (include_viscosity) {
    setPermissions(readIfSet("species:e:density", Regions::Interior));
  }
  // Add phi to restart files so that the value in the boundaries
  // is restored on restart. This is done even when phi is not evolving,
  // so that phi can be saved and re-loaded

  // Coefficients for Boussinesq solve
  boussinesq_coeff = average_atomic_mass / SQ(coord->Bxy());
  reynolds_coefficients = bout::ion_energy::reynoldsMetricCoefficients(boussinesq_coeff);

  phiSolver = Laplacian::create(&options["laplacian"]);
  phiSolver->setCoefC(boussinesq_coeff);

  if (radial_model.walls) {
    phiSolver->setInnerBoundaryFlags(radial_model.phi_flags);
    phiSolver->setOuterBoundaryFlags(radial_model.phi_flags);
  }

  // Coupled pressure/potential recovery settings
  phi_max_iterations = options["phi_max_iterations"]
                           .doc("Maximum number of finite-volume potential iterations")
                           .withDefault(100);
  phi_atol = options["phi_atol"]
                 .doc("Potential scale in the relative recovery criterion [eV]")
                 .withDefault(1.0)
             / Tnorm;
  phi_rtol =
      options["phi_rtol"].doc("phi convergence relative tolerance").withDefault(1e-3);
  vort_atol = options["vort_atol"]
                  .doc("Absolute tolerance on the normalised forward vorticity residual")
                  .withDefault(1e-10);

  // Read curvature vector
  try {
    Curlb_B.covariant = false; // Contravariant
    mesh->get(Curlb_B, "bxcv");

  } catch (BoutException& e) {
    try {
      // May be 2D, reading as 3D
      Vector2D curv2d;
      curv2d.covariant = false;
      mesh->get(curv2d, "bxcv");
      Curlb_B = curv2d;
    } catch (BoutException& e) {
      output_warn.write("No curvature vector in input grid");
      Curlb_B = 0.0;
    }
  }

  if (Options::root()["mesh"]["paralleltransform"]["type"].as<std::string>()
      == "shifted") {
    Field2D I;
    mesh->get(I, "sinty");
    Curlb_B.z += I * Curlb_B.x;
  }

  Curlb_B.x /= Bnorm;
  Curlb_B.y *= SQ(Lnorm);
  Curlb_B.z *= SQ(Lnorm);

  Curlb_B *= 2. / coord->Bxy();

  diagnose =
      options["diagnose"].doc("Output additional diagnostics?").withDefault<bool>(false);
}

void EvolveVorticityIonEnergy::acceptWarmStart() {
  if (!have_solution) {
    return; // No converged recovery yet: keep the initial or restart guess
  }
  phi_guess = copy(phi);
  ion_KE_perp_guess = copy(ion_KE_perp);
}

void EvolveVorticityIonEnergy::transform_impl(GuardedOptions& state) {
  mesh->communicate(Vort, Eion);
  vorticity = copy(Vort);

  // Start from the guess stored at the last output state, not from whatever
  // the previous RHS call (possibly a Jacobian probe or rejected solver evaluation) left.
  // This keeps the RHS a deterministic function of the state within a step.
  phi = copy(phi_guess);

  if (radial_model.walls) {
    // Assignment preserves vorticity's independently owned boundary rules.
    vorticity.applyBoundary();
    bout::ion_energy::applyVorticityWallBoundary(vorticity,
                                                 radial_model.vorticity_wall_zero);
    mesh->communicate(vorticity);
    applyPotentialBoundary();
  }

  auto allspecies = state["species"];
  auto electrons = allspecies["e"];

  // From vorticity and ion energy E we will recover
  // the potential phi (ϕ) and ion pressure p_i
  //
  // With Boussinesq approximation this can be solved once for Φ but
  // applying boundary conditions on pressure and phi separately
  // requires iteration.

  // Calculate K_||i = ( m_i n v_||i )^2 / (2 m_i n)

  ion_KE_parallel = 0.0;
  Field3D Ni_sum = 0.0;
  Field3D mass_charge_density = 0.0;
  Field3D mass_density = 0.0;
  forAllIons(allspecies, [&](GuardedOptions species) {
    const BoutReal AA = get<BoutReal>(species["AA"]);
    const BoutReal Z = get<BoutReal>(species["charge"]);
    auto N = GET_NOBOUNDARY(Field3D, species["density"]);
    Ni_sum += N;
    mass_charge_density += (AA / Z) * N;
    mass_density += AA * N;

    if (!species.isSet("velocity")) {
      return;
    }
    auto V = GET_NOBOUNDARY(Field3D, species["velocity"]);

    ion_KE_parallel += 0.5 * AA * N * SQ(V);
  });

  const Field3D Ni_sum_floor = softFloor(Ni_sum, density_floor);
  // With a common temperature, psi = Ti * sum_s(A_s n_s / Z_s) / Abar.
  // Use the same density floor as the temperature recovery.
  const Field3D pressure_weight =
      mass_charge_density / (average_atomic_mass * Ni_sum_floor);
  // Initial pressure: this state's E and K_|| with K⊥ from the last output
  // state. Equivalent to reusing the previous Pi_sum corrected for the change in
  // (E - K_||), so only the change in K⊥ is left for the iteration. Before the
  // first output the guess is K⊥ = 0.
  // Pressure guards must be consistent before the first Phi correction.
  Pi_sum = floor((2. / 3) * (Eion - ion_KE_parallel - ion_KE_perp_guess), 0.0);
  bout::ion_energy::applyWallPressureBoundary(Pi_sum);
  psi = pressure_weight * Pi_sum;
  bout::ion_energy::applyWallPressureBoundary(psi);
  const BoutReal vort_tolerance = vort_atol + (phi_rtol * max(abs(Vort), true));

  int iteration = 0;
  BoutReal phi_err;
  do {
    // Both potentials must converge: phi can be unchanged on the first
    // iteration even though the recovered pressure potential has changed.
    const Field3D phi_last = copy(phi);
    const Field3D psi_last = copy(psi);

    // Solve
    // ∇·(c∇⊥Φ) = ω for Φ = ϕ + psi
    // with c = rho_0 / B^2 and psi = sum_s(A_s p_s / (Abar Z_s))

    mesh->communicate(phi, psi);
    applyPotentialBoundary();
    // The spectral solver is only an approximate inverse. Defect correction
    // converges to the same finite-volume operator used for kinetic energy.
    phi = correctPotential(Vort, phi + psi) - psi;
    applyPotentialBoundary();

    // With the Boussinesq approximation the perpendicular flow velocity is
    // v⊥0 = b × ∇Φ / B = v_E + v∗
    // The mass density is constant so kinetic energy K⊥ = 1/2 c |∇⊥Φ|^2,

    ion_KE_perp = calculateKEperp(phi, psi);

    // Calculate pressure
    // E = Cv * P + (1/2) m n v^2 + K_\perp
    // Ensure that P does not go negative.
    Pi_sum = floor((2. / 3) * (Eion - ion_KE_parallel - ion_KE_perp), 0.0);

    bout::ion_energy::applyWallPressureBoundary(Pi_sum);
    psi = pressure_weight * Pi_sum;
    bout::ion_energy::applyWallPressureBoundary(psi);
    mesh->communicate(phi, psi);
    applyPotentialBoundary();

    phi_err = max((abs(phi - phi_last)) / (abs(phi) + phi_atol), true);
    const BoutReal psi_err = max(abs(psi - psi_last) / (abs(psi) + phi_atol), true);
    phi_err = std::max(phi_err, psi_err);
    const BoutReal vort_err =
        max(abs(Vort - bout::ion_energy::Div_a_Grad_perp_XZ(boussinesq_coeff, phi + psi)),
            true);
    ++iteration;

    if (phi_err <= phi_rtol && vort_err <= vort_tolerance) {
      break;
    }
    if (iteration >= phi_max_iterations) {
      throw BoutException("Potential recovery did not converge in {} iterations: "
                          "potential error {}, vorticity residual {} (tolerance {})",
                          iteration, phi_err, vort_err, vort_tolerance);
    }
  } while (true);
  phi_iterations = iteration;
  mesh->communicate(phi, Pi_sum, psi);
  applyPotentialBoundary();

  // Free parallel boundaries, as in the vorticity component. The spectral
  // inverse does not supply these, but electron compression uses Grad(phi).
  if (mesh->ystart > 0) {
    Field3D phi_fa = toFieldAligned(phi);
    const int lower_neighbour = std::min(mesh->ystart + 1, mesh->yend);
    const int upper_neighbour = std::max(mesh->yend - 1, mesh->ystart);
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        phi_fa(r.ind, mesh->ystart - 1, z) =
            (2 * phi_fa(r.ind, mesh->ystart, z)) - phi_fa(r.ind, lower_neighbour, z);
      }
    }
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        phi_fa(r.ind, mesh->yend + 1, z) =
            (2 * phi_fa(r.ind, mesh->yend, z)) - phi_fa(r.ind, upper_neighbour, z);
      }
    }
    phi = fromFieldAligned(phi_fa);
  }
  applyPotentialBoundary();
  have_solution = true; // phi and ion_KE_perp are converged for this state
  if (radial_model.walls) {
    computeWallRates();
  }

  auto fields = state["fields"];
  set(fields["vorticity"], vorticity);
  set(fields["phi"], phi);

  //////////////////////////////////////////////////////////
  // Calculate common ion temperature
  // Set temperature and pressure of each ion species

  Ti = Pi_sum / Ni_sum_floor;
  forAllIons(allspecies, [&](GuardedOptions species) {
    set(species["temperature"], Ti);
    auto Ni = GET_NOBOUNDARY(Field3D, species["density"]);
    set(species["pressure"], Ni * Ti);
  });

  //////////////////////////////////////////////////////////
  // Reynolds stress

  const Field3D Phi = phi + psi;
  const auto reynolds = bout::ion_energy::reynoldsFaces(
      Field3D{boussinesq_coeff}, vorticity, phi, psi, reynolds_coefficients);
  DivJReynolds = -bout::ion_energy::divergence(reynolds);
  if (poloidal_flows) {
    DivJReynolds -= Div_n_bxGrad_f_B_XPPM(vorticity, phi, false, true)
                    - Div_n_bxGrad_f_B_XPPM(vorticity, phi, false, false);
  }
  if (diagnose) {
    reynolds_work = bout::ion_energy::faceWork(reynolds, Phi);
  }

  //////////////////////////////////////////////////////////
  // Diamagnetic drift

  // Electron diamagnetic current
  const auto Pe = GET_NOBOUNDARY(Field3D, electrons["pressure"]);
  const Vector3D Jdia_e = Pe * Curlb_B;

  // Ion diamagnetic current
  const Vector3D Jdia_i = Pi_sum * Curlb_B;

  DivJdia = Div(Jdia_e + Jdia_i);
  const Field3D DivJdia_e = Div(Jdia_e);

  // Diamagnetic electrical terms belong to recovery. Parallel electrical
  // terms require boundary-conditioned momentum and potential in finally().
  DivEdiamagnetic = Div(phi * Jdia_i);
  DivEelectric = DivEdiamagnetic;

  //////////////////////////////////////////////////////////
  // ExB drift

  // ExB drift of total ion energy. This includes the Reynolds work that
  // is the advection of K_perp by ExB for the Reynolds current used here.
  DivEExB = Div_n_bxGrad_f_B_XPPM(Eion, phi, false, poloidal_flows, true);

  // The total electron energy exchange is phi * Div(J_e). evolve_pressure
  // evolves thermal energy, so only ExB compression belongs in energy_source.
  // Parallel electrical work is already in the electron momentum equation.
  phi_Div_Jdia_e = phi * DivJdia_e;
  phi_Div_Je = phi_Div_Jdia_e;
  electron_exb_work = -(Jdia_e * Grad(phi));
  add(electrons["energy_source"], electron_exb_work);

  //////////////////////////////////////////////////////////
  // Kinematic viscosity

  // The polarisation energy transport already carries viscous drift work and
  // enthalpy. Do not add a separate Phi-weighted viscous energy flux: that would
  // count the work twice. Ion heating follows from total-energy recovery;
  // the discrete stress-work diagnostic supplies the electron heating share.
  if (include_viscosity) {
    calculateViscosity();
    // Match vorticity's mass-density partition, including electrons when
    // their mass and density are available. Ion heating is already internal
    // to Eion; only transfer the electron share out of the ion reservoir.
    if (electrons.isSet("AA") && electrons.isSet("density")) {
      const Field3D electron_mass_density =
          get<BoutReal>(electrons["AA"]) * GET_NOBOUNDARY(Field3D, electrons["density"]);
      electron_viscous_heating = viscous_heating * electron_mass_density
                                 / floor(mass_density + electron_mass_density,
                                         average_atomic_mass * density_floor);
      add(electrons["energy_source"], electron_viscous_heating);
    }
  }

  // Physics inputs for conservative drift, independent of diagnostics.
  set(fields["ion_boussinesq_coefficient"], boussinesq_coeff);
  set(fields["ion_pressure_potential"], psi);
  set(fields["ion_mass_density"], mass_density);
  set(fields["ion_reynolds_x"], reynolds.x);
  set(fields["ion_reynolds_z"], reynolds.z);
  set(fields["DivJReynolds"], DivJReynolds);
  set(fields["DivJdia"], DivJdia);
  Field3D wall_rate = zeroFrom(phi);
  writeRadialValues(wall_rate, wall_rate_inner, wall_rate_outer);
  set(fields["ion_wall_rate"], wall_rate);
}

void EvolveVorticityIonEnergy::finally(const Options& state) {
  const Field3D DivEpol =
      get<Field3D>(state["fields"]["ion_polarisation_energy_divergence"]);
  const Field3D current_rhs = get<Field3D>(state["fields"]["ion_current_rhs"]);

  // The completed potential includes the sheath's guards. Keep it for output
  // too, so saved target work can be checked against the numerical expression.
  phi = get<Field3D>(state["fields"]["phi"]);
  ion_phi_dissipation = zeroFrom(phi);
  if (phi_dissipation) {
    if (!state.isSet("sound_speed")) {
      throw BoutException("evolve_vorticity_ion_energy:phi_dissipation requires "
                          "completed sound_speed; enable the sound_speed component");
    }
    // Same MC face flux and default physical-Y flux treatment as Vorticity.
    // Sound speed is constructed after recovery publishes ion pressure.
    const Field3D damping = -FV::Div_par(-phi, 0.0, get<Field3D>(state["sound_speed"]));
    BOUT_FOR(i, Vort.getRegion("RGN_NOBNDRY")) { ion_phi_dissipation[i] = damping[i]; }
  }
  Jpar_e = zeroFrom(phi);
  Jpar_i = zeroFrom(phi);
  const auto& electrons = state["species"]["e"];
  if (electrons.isSet("momentum")) {
    Jpar_e = (get<BoutReal>(electrons["charge"]) / get<BoutReal>(electrons["AA"]))
             * get<Field3D>(electrons["momentum"]);
  }
  forAllIons(state["species"], [&](const Options& species) {
    if (species.isSet("momentum")) {
      Jpar_i += (get<BoutReal>(species["charge"]) / get<BoutReal>(species["AA"]))
                * get<Field3D>(species["momentum"]);
    }
  });
  mesh->communicate(Jpar_e, Jpar_i);
  DivEelectric = DivEdiamagnetic + Div_par(phi * Jpar_i);
  phi_Div_Je = phi_Div_Jdia_e + phi * Div_par(Jpar_e);

  //////////////////////////////////////////////////////////
  // Vorticity equation

  if (radial_model.walls) {
    ddt(Vort) = 0.0;
    BOUT_FOR(i, Vort.getRegion("RGN_NOBNDRY")) {
      ddt(Vort)[i] = current_rhs[i] + DivJviscous[i] + ion_phi_dissipation[i];
    }
    writeWallRates();
  } else {
    ddt(Vort) = current_rhs + DivJviscous;
    BOUT_FOR(i, Vort.getRegion("RGN_NOBNDRY")) { ddt(Vort)[i] += ion_phi_dissipation[i]; }
  }

  //////////////////////////////////////////////////////////
  // Ion energy equation

  // Each species carries its own parallel enthalpy and kinetic energy.
  // Read these from the completed state so boundary values are available.
  DivEparallel = 0.0;
  forAllIons(state["species"], [&](const Options& species) {
    if (!species.isSet("velocity")) {
      return;
    }
    const Field3D N = get<Field3D>(species["density"]);
    const Field3D V = get<Field3D>(species["velocity"]);
    const Field3D P = get<Field3D>(species["pressure"]);
    const BoutReal AA = get<BoutReal>(species["AA"]);
    const Field3D energy = (5. / 2) * P + 0.5 * AA * N * SQ(V);
    const Field3D fastest_wave = state.isSet("fastest_wave")
                                     ? get<Field3D>(state["fastest_wave"])
                                     : sqrt(get<Field3D>(species["temperature"]) / AA);
    Field3D flow_ylow;
    DivEparallel += FV::Div_par_mod<hermes::Limiter>(energy, V, fastest_wave, flow_ylow);
  });

  Eion_src = sumIonEnergySources(state["species"]);

  ddt(Eion) = -(DivEExB + DivEelectric + DivEparallel + DivEpol) - phi_Div_Je
              - electron_viscous_heating + Eion_src;
}

void EvolveVorticityIonEnergy::calculateViscosity() {
  const Field3D nu{viscosity};
  const Field3D Phi = phi + psi;
  auto flux = bout::ion_energy::perpendicularFaces(nu, vorticity);
  // The stress diagnostic has the same wall closure as the vorticity rate.
  // Its Phi-weighted flux is diagnostic only; DivEpol carries the drift work
  // for both open and closed viscous wall faces.
  if (!radial_model.viscous_wall_flux) {
    bout::ion_energy::closeRadialFaces(flux);
  }
  DivJviscous = bout::ion_energy::divergence(flux);
  DivEviscous = zeroFrom(vorticity);
  electron_viscous_heating = 0.0;
  const auto* coord = mesh->getCoordinates();
  BOUT_FOR(i, vorticity.getRegion("RGN_NOBNDRY")) {
    DivEviscous[i] =
        0.5
        * ((Phi[i] + Phi[i.xp()]) * flux.x[i] - (Phi[i] + Phi[i.xm()]) * flux.x[i.xm()]
           + (Phi[i] + Phi[i.zp()]) * flux.z[i] - (Phi[i] + Phi[i.zm()]) * flux.z[i.zm()])
        / (coord->J()[i] * coord->dx()[i] * coord->dy()[i] * coord->dz()[i]);
  }
  viscous_heating = Phi * DivJviscous - DivEviscous;
}

void EvolveVorticityIonEnergy::outputVars(Options& state) {
  // Outputs follow an RHS evaluation at the accepted output state.
  // Update the recovery guess here so it remains fixed between outputs.
  acceptWarmStart();

  set_with_attrs(
      state["evolve_vorticity_ion_energy_version"], 7,
      {{"description", "Component equation and discretisation version: "
                       "MC full-vorticity advection plus compact central R; "
                       "Phi-weighted polarisation work with wall-uniform transport; "
                       "viscous drift energy carried by polarisation transport; "
                       "conservative polarisation owned by polarisation_drift; "
                       "boundary-conditioned parallel/extra/collisional currents; "
                       "completed-state parallel electrical work; "
                       "optional collective-sound-speed potential damping; "
                       "independent of restart layout"}});
  set_with_attrs(
      state["ion_phi_dissipation_enabled"], phi_dissipation,
      {{"description", "Resolved numerical parallel potential damping option"}});

  // Normalisations
  auto Nnorm = get<BoutReal>(state["Nnorm"]);
  auto Tnorm = get<BoutReal>(state["Tnorm"]);
  auto Omega_ci = get<BoutReal>(state["Omega_ci"]);
  auto rho_s0 = get<BoutReal>(state["rho_s0"]);

  BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation

  state["Eion"].setAttributes({{"time_dimension", "t"},
                               {"units", "J/m^3"},
                               {"conversion", Pnorm},
                               {"standard_name", "E_ion"},
                               {"long_name", "Ion energy density"},
                               {"species", "all ions"},
                               {"source", "evolve_vorticity_ion_energy"}});

  set_with_attrs(state["Tion"], Ti,
                 {{"time_dimension", "t"},
                  {"units", "eV"},
                  {"conversion", Tnorm},
                  {"standard_name", "temperature"},
                  {"long_name", "Common ion temperature"},
                  {"species", "all ions"},
                  {"source", "evolve_vorticity_ion_energy"}});

  if (radial_model.walls) {
    set(state["Vort"], vorticity);
  }
  state["Vort"].setAttributes({{"time_dimension", "t"},
                               {"units", "C m^-3"},
                               {"conversion", SI::qe * Nnorm},
                               {"long_name", "vorticity"},
                               {"source", "evolve_vorticity_ion_energy"}});

  set_with_attrs(state["phi"], phi,
                 {{"time_dimension", "t"},
                  {"units", "V"},
                  {"conversion", Tnorm},
                  {"standard_name", "potential"},
                  {"long_name", "plasma potential"},
                  {"source", "evolve_vorticity_ion_energy"}});

  if (diagnose) {
    set_with_attrs(state["ion_phi_dissipation"], ion_phi_dissipation,
                   {{"time_dimension", "t"},
                    {"units", "A m^-3"},
                    {"conversion", SI::qe * Nnorm * Omega_ci},
                    {"long_name", "Numerical parallel potential vorticity rate"},
                    {"description", "Minus FV::Div_par(-phi,0,sound_speed); excluded "
                                    "from chi and physical currents; no explicit "
                                    "heating or external power sink"}});
    set_with_attrs(state["DivJReynolds"], DivJReynolds,
                   {{"time_dimension", "t"},
                    {"units", "A m^-3"},
                    {"conversion", SI::qe * Nnorm * Omega_ci},
                    {"long_name", "Reynolds current divergence"}});
    set_with_attrs(state["ddt(Vort)"], ddt(Vort),
                   {{"time_dimension", "t"},
                    {"units", "A m^-3"},
                    {"conversion", SI::qe * Nnorm * Omega_ci},
                    {"long_name", "Completed vorticity rate in physical cells"}});
    for (const auto& [name, current] :
         {std::pair{"ion_parallel_current", &Jpar_i},
          std::pair{"electron_parallel_current", &Jpar_e}}) {
      set_with_attrs(state[name], *current,
                     {{"time_dimension", "t"},
                      {"units", "A m^-2"},
                      {"conversion", SI::qe * Nnorm * rho_s0 * Omega_ci},
                      {"source", "evolve_vorticity_ion_energy"},
                      {"long_name", "Completed parallel current with sheath guards"}});
    }
    if (radial_model.walls) {
      Field3D wall_rate = zeroFrom(phi);
      writeRadialValues(wall_rate, wall_rate_inner, wall_rate_outer);
      set_with_attrs(state["phi_wall_rate"], wall_rate,
                     {{"time_dimension", "t"},
                      {"units", "V s^-1"},
                      {"conversion", Tnorm * Omega_ci},
                      {"long_name", "Radial midpoint potential rate"},
                      {"description", "Wall rates in physical radial guard cells; "
                                      "zero in interior and Y boundary cells"}});
    }
    if (include_viscosity) {
      set_with_attrs(state["kinematic_viscosity"], viscosity,
                     {{"units", "m^2 s^-1"},
                      {"conversion", SQ(rho_s0) * Omega_ci},
                      {"long_name", "Kinematic viscosity"}});
      set_with_attrs(state["kinematic_viscous_heating"], viscous_heating,
                     {{"time_dimension", "t"},
                      {"units", "W m^-3"},
                      {"conversion", Pnorm * Omega_ci},
                      {"long_name", "Discrete kinematic viscous heating"}});
      set_with_attrs(state["electron_viscous_heating"], electron_viscous_heating,
                     {{"time_dimension", "t"},
                      {"units", "W m^-3"},
                      {"conversion", Pnorm * Omega_ci},
                      {"long_name", "Electron share of kinematic viscous heating"}});
      set_with_attrs(state["DivJviscous"], DivJviscous,
                     {{"time_dimension", "t"},
                      {"units", "A m^-3"},
                      {"conversion", SI::qe * Nnorm * Omega_ci},
                      {"long_name", "Kinematic viscous current divergence"}});
      set_with_attrs(
          state["DivEviscous"], DivEviscous,
          {{"time_dimension", "t"},
           {"units", "W m^-3"},
           {"conversion", Pnorm * Omega_ci},
           {"long_name", "Diagnostic Phi-weighted viscous stress flux divergence"},
           {"description", "Viscous drift energy is carried by "
                           "polarisation transport"}});
    }
    set_with_attrs(state["ion_pressure_potential"], psi,
                   {{"time_dimension", "t"}, {"units", "V"}, {"conversion", Tnorm}});
    set_with_attrs(state["ion_phi_iterations"], phi_iterations,
                   {{"time_dimension", "t"},
                    {"description", "Iterations used in the last coupled pressure and "
                                    "potential recovery"}});
    set_with_attrs(state["ion_KE_parallel"], ion_KE_parallel,
                   {{"time_dimension", "t"}, {"units", "J m^-3"}, {"conversion", Pnorm}});
    set_with_attrs(state["ion_KE_perp"], ion_KE_perp,
                   {{"time_dimension", "t"}, {"units", "J m^-3"}, {"conversion", Pnorm}});
    set_with_attrs(
        state["ddt(Eion)"], ddt(Eion),
        {{"time_dimension", "t"}, {"units", "W m^-3"}, {"conversion", Pnorm * Omega_ci}});
    set_with_attrs(
        state["DivEExB"], DivEExB,
        {{"time_dimension", "t"}, {"units", "W m^-3"}, {"conversion", Pnorm * Omega_ci}});
    set_with_attrs(
        state["reynolds_work"], reynolds_work,
        {{"time_dimension", "t"}, {"units", "W m^-3"}, {"conversion", Pnorm * Omega_ci}});
    set_with_attrs(
        state["DivEelectric"], DivEelectric,
        {{"time_dimension", "t"}, {"units", "W m^-3"}, {"conversion", Pnorm * Omega_ci}});
    set_with_attrs(
        state["DivEparallel"], DivEparallel,
        {{"time_dimension", "t"}, {"units", "W m^-3"}, {"conversion", Pnorm * Omega_ci}});
    set_with_attrs(
        state["phi_Div_Je"], phi_Div_Je,
        {{"time_dimension", "t"}, {"units", "W m^-3"}, {"conversion", Pnorm * Omega_ci}});
    set_with_attrs(
        state["electron_exb_work"], electron_exb_work,
        {{"time_dimension", "t"}, {"units", "W m^-3"}, {"conversion", Pnorm * Omega_ci}});
    set_with_attrs(
        state["Eion_src"], Eion_src,
        {{"time_dimension", "t"}, {"units", "W m^-3"}, {"conversion", Pnorm * Omega_ci}});
  }
}

void EvolveVorticityIonEnergy::initialiseBoundaryState() {
  if (!radial_model.walls) {
    return;
  }
  BOUT_FOR(i, Vort.getRegion("RGN_BNDRY")) { Vort[i] = 0.0; }
  const auto initialise = [&](int interior, int guard, int begin, int end,
                              bool average_y) {
    std::vector<BoutReal> midpoints(static_cast<size_t>(mesh->LocalNy), 0.0);
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      midpoints[static_cast<size_t>(y)] =
          0.5 * (meanZ(phi, guard, y) + meanZ(phi, interior, y));
    }
    if (average_y) {
      const BoutReal average = fluxSurfaceAverage(midpoints, interior, guard);
      std::fill(midpoints.begin(), midpoints.end(), average);
    }
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      for (int x = begin; x < end; ++x) {
        for (int z = 0; z < mesh->LocalNz; ++z) {
          Vort(x, y, z) = midpoints[static_cast<size_t>(y)];
        }
      }
    }
  };
  if (mesh->firstX()) {
    initialise(mesh->xstart, mesh->xstart - 1, 0, mesh->xstart,
               radial_model.core_flux_surface_average);
  }
  if (mesh->lastX()) {
    initialise(mesh->xend, mesh->xend + 1, mesh->xend + 1, mesh->LocalNx, false);
  }
}

void EvolveVorticityIonEnergy::applyPotentialBoundary() {
  if (!radial_model.walls) {
    return;
  }
  const auto apply = [&](int interior, int guard, int begin, int end, bool average_y) {
    std::vector<BoutReal> midpoints(static_cast<size_t>(mesh->LocalNy), 0.0);
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      midpoints[static_cast<size_t>(y)] = meanZ(Vort, guard, y);
    }
    if (average_y) {
      const BoutReal average = fluxSurfaceAverage(midpoints, interior, guard);
      std::fill(midpoints.begin(), midpoints.end(), average);
    }
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        phi(guard, y, z) = 2.0 * midpoints[static_cast<size_t>(y)] - phi(interior, y, z);
        for (int x = begin; x < end; ++x) {
          phi(x, y, z) = phi(guard, y, z);
        }
      }
    }
  };
  if (mesh->firstX()) {
    apply(mesh->xstart, mesh->xstart - 1, 0, mesh->xstart - 1,
          radial_model.core_flux_surface_average);
  }
  if (mesh->lastX()) {
    apply(mesh->xend, mesh->xend + 1, mesh->xend + 2, mesh->LocalNx, false);
  }
}

void EvolveVorticityIonEnergy::computeWallRates() {
  // Reset every RHS evaluation; no rate survives from a previous solver state.
  std::fill(wall_rate_inner.begin(), wall_rate_inner.end(), 0.0);
  std::fill(wall_rate_outer.begin(), wall_rate_outer.end(), 0.0);
  if (!radial_model.walls || !radial_model.evolve_wall_potential) {
    return;
  }
  const auto* coord = mesh->getCoordinates();
  const auto compute = [&](int interior, int guard, bool average_y,
                           std::vector<BoutReal>& rates) {
    BoutReal core_midpoint = 0.0;
    if (average_y) {
      std::vector<BoutReal> midpoints(static_cast<size_t>(mesh->LocalNy), 0.0);
      for (int y = mesh->ystart; y <= mesh->yend; ++y) {
        midpoints[static_cast<size_t>(y)] = meanZ(Vort, guard, y);
      }
      core_midpoint = fluxSurfaceAverage(midpoints, interior, guard);
    }
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      const BoutReal midpoint = average_y ? core_midpoint : meanZ(Vort, guard, y);
      const BoutReal derivative = bout::ion_energy::wallInwardDerivative(
          midpoint, meanZ(phi, interior, y), meanZ(coord->dx(), interior, y));
      // For a 2D metric the mean g13 contribution vanishes on periodic Z.
      // With BOUT_USE_METRIC_3D we use Z-averaged geometry, averaging g11
      // before taking the square root to estimate the DC normal derivative.
      const BoutReal normal_metric =
          sqrt(0.5 * (meanZ(coord->g11(), guard, y) + meanZ(coord->g11(), interior, y)));
      rates[static_cast<size_t>(y)] =
          radial_model.phi_boundary_velocity * normal_metric * derivative;
    }
    if (average_y) {
      const BoutReal average = fluxSurfaceAverage(rates, interior, guard);
      // Y guard entries stay zero, including the DC values passed to chi.
      std::fill(rates.begin() + mesh->ystart, rates.begin() + mesh->yend + 1, average);
    }
  };
  if (mesh->firstX()) {
    compute(mesh->xstart, mesh->xstart - 1, radial_model.core_flux_surface_average,
            wall_rate_inner);
  }
  if (mesh->lastX()) {
    compute(mesh->xend, mesh->xend + 1, false, wall_rate_outer);
  }
}

void EvolveVorticityIonEnergy::writeWallRates() {
  writeRadialValues(ddt(Vort), wall_rate_inner, wall_rate_outer);
}

void EvolveVorticityIonEnergy::restartVars(Options& state) {
  // The first call is the postInit hook, after solver variables were loaded.
  if (restore_phi_on_restart) {
    restore_phi_on_restart = false;
    if (state.isSet("phi")) {
      phi = state["phi"].as<Field3D>();
      phi_guess = copy(phi);
    }
    if (radial_model.walls) {
      const int version = state["ion_energy_boundary_state_version"].withDefault(0);
      if (version < 0 || version > 1) {
        throw BoutException("Unsupported ion energy boundary state version {}", version);
      }
      if (version == 0) {
        initialiseBoundaryState();
      }
    }
  }
  if (radial_model.walls) {
    set_with_attrs(state["ion_energy_boundary_state_version"], 1,
                   {{"description", "Restart Vort physical radial boundary cells store "
                                    "phi midpoints in normalised potential units"}});
  }
  set_with_attrs(
      state["phi"], phi,
      {{"long_name", "plasma potential"}, {"source", "evolve_vorticity_ion_energy"}});
}

Field3D EvolveVorticityIonEnergy::sumIonEnergySources(const Options& allspecies) const {
  Field3D src = 0.0;
  forAllIons(allspecies, [&](const Options& species) {
    if (species.isSet("energy_source")) {
      src += get<Field3D>(species["energy_source"]); // Thermal energy only
    }
    if (species.isSet("velocity")) {
      const Field3D V = get<Field3D>(species["velocity"]);
      if (species.isSet("momentum_source")) {
        src += V * get<Field3D>(species["momentum_source"]);
      }
      if (species.isSet("density_source")) {
        const BoutReal AA = get<BoutReal>(species["AA"]);
        src -= 0.5 * AA * SQ(V) * get<Field3D>(species["density_source"]);
      }
      // Polarisation advection is already present in the total-energy flux.
      // Remove its local kinetic-energy conversion using the completed V.

      if (species.isSet("polarisation_momentum_source")) {
        src -= V * get<Field3D>(species["polarisation_momentum_source"]);
      }
      if (species.isSet("polarisation_density_source")) {
        const BoutReal AA = get<BoutReal>(species["AA"]);
        src += 0.5 * AA * SQ(V) * get<Field3D>(species["polarisation_density_source"]);
      }
    }
  });
  return src;
}

// Use exactly the same face fluxes as the forward vorticity operator.
Field3D EvolveVorticityIonEnergy::calculateKEperp(const Field3D& phi,
                                                  const Field3D& psi) const {
  return bout::ion_energy::kineticEnergyPerp(boussinesq_coeff, phi + psi);
}

Field3D EvolveVorticityIonEnergy::correctPotential(const Field3D& rhs,
                                                   const Field3D& potential) const {
  return correctPotential(rhs, potential, *phiSolver);
}

Field3D EvolveVorticityIonEnergy::correctPotential(const Field3D& rhs,
                                                   const Field3D& potential,
                                                   Laplacian& approximate_inverse) const {
  // Homogeneous correction boundaries preserve the separately imposed phi
  // boundaries, including INVERT_SET values supplied by a restart.
  Field3D residual{zeroFrom(rhs)};
  const Field3D forward =
      bout::ion_energy::Div_a_Grad_perp_XZ(boussinesq_coeff, potential);
  BOUT_FOR(i, residual.getRegion("RGN_NOBNDRY")) {
    residual[i] = (rhs[i] - forward[i]) / boussinesq_coeff[i];
  }
  const Field3D approximate = approximate_inverse.solve(residual, zeroFrom(potential));
  // Spectral solvers can leave Y guard cells uninitialised. Only copy the
  // X-Z planes they solve; communication supplies inter-process Y guards.
  Field3D correction{zeroFrom(potential)};
  BOUT_FOR(i, correction.getRegion("RGN_NOY")) { correction[i] = approximate[i]; }
  mesh->communicate(correction);
  return potential + correction;
}
