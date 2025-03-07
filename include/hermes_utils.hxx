#pragma once
#ifndef HERMES_UTILS_H
#define HERMES_UTILS_H

#include <bout/mesh.hxx>

inline BoutReal floor(BoutReal value, BoutReal min) {
  if (value < min)
    return min;
  return value;
}

template<typename T, typename = bout::utils::EnableIfField<T>>
inline T clamp(const T& var, BoutReal lo, BoutReal hi, const std::string& rgn = "RGN_ALL") {
  checkData(var);
  T result = copy(var);

  BOUT_FOR(d, var.getRegion(rgn)) {
    if (result[d] < lo) {
      result[d] = lo;
    } else if (result[d] > hi) {
      result[d] = hi;
    }
  }

  return result;
}

template<typename T, typename = bout::utils::EnableIfField<T>>
Ind3D indexAt(const T& f, int x, int y, int z) {
  int ny = f.getNy();
  int nz = f.getNz();
  return Ind3D{(x * ny + y) * nz + z, ny, nz};
}

/// Utility for iterating over boundary points
///
/// The function should take 3 arguments:
///   void func (index_bndry, index_domain, index_domain2, sign)
///
/// index_bndry is the Ind3D index of a point in the boundary
/// index_domain is the neighboring point in the domain
/// index_domain2 is the next point into the domain
/// sign is -1 for lower Y boundaries, +1 for upper Y.
template <typename Function>
void iterateBoundaries(Function func) {
  const int ny = bout::globals::mesh->LocalNy;
  const int nz = bout::globals::mesh->LocalNz;

  // Lower boundary
  const int ystart = bout::globals::mesh->ystart;
  for (RangeIterator r = bout::globals::mesh->iterateBndryLowerY(); !r.isDone(); r++) {
    for (int jz = 0; jz < nz; jz++) {
      const auto i = Ind3D{(r.ind * ny + ystart) * nz + jz, ny, nz};
      const auto im = i.ym(); // In boundary
      const auto ip = i.yp(); // Away from boundary
      func(im, i, ip, -1);
    }
  }

  // Upper boundary
  const int yend = bout::globals::mesh->yend;
  for (RangeIterator r = bout::globals::mesh->iterateBndryUpperY(); !r.isDone(); r++) {
    for (int jz = 0; jz < nz; jz++) {
      auto i = Ind3D{(r.ind * ny + yend) * nz + jz, ny, nz};
      auto im = i.ym();
      auto ip = i.yp(); // Into boundary
      func(ip, i, im, +1);
    }
  }
}

#endif // HERMES_UTILS_H
