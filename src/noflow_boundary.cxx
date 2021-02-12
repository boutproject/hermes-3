#include "../include/noflow_boundary.hxx"

#include "bout/mesh.hxx"
using bout::globals::mesh;

namespace {
Ind3D indexAt(const Field3D& f, int x, int y, int z) {
  int ny = f.getNy();
  int nz = f.getNz();
  return Ind3D{(x * ny + y) * nz + z, ny, nz};
}
}

void NoFlowBoundary::transform(Options& state) {
  AUTO_TRACE();

  // Make sure that the state has been set for this species
  ASSERT1(state.isSection(name));

  Options& species = state[name];

  // Apply zero-gradient boundary conditions to state variables
  for (const std::string field : {"density", "temperature", "pressure"}) {
    if (!species.isSet(field)) {
      continue; // Skip variables which are not set
    }

    Field3D var = get<Field3D>(species[field]);

    if (noflow_lower_y) {
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          auto i = indexAt(var, r.ind, mesh->ystart, jz);
          auto im = i.ym();

          var[im] = var[i];
        }
      }
    }

    if (noflow_upper_y) {
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          auto i = indexAt(var, r.ind, mesh->yend, jz);
          auto ip = i.yp();

          var[ip] = var[i];
        }
      }
    }

    set<Field3D>(species[field], var);
  }

  // Apply zero-value boundary conditions to flows
  for (const std::string field : {"velocity", "momentum"}) {
    if (!species.isSet(field)) {
      continue; // Skip variables which are not set
    }

    Field3D var = get<Field3D>(species[field]);

    if (noflow_lower_y) {
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          auto i = indexAt(var, r.ind, mesh->ystart, jz);
          auto im = i.ym();

          var[im] = - var[i];
        }
      }
    }

    if (noflow_upper_y) {
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          auto i = indexAt(var, r.ind, mesh->yend, jz);
          auto ip = i.yp();

          var[ip] = - var[i];
        }
      }
    }

    set<Field3D>(species[field], var);
  }
}
