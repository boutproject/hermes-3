#pragma once
#ifndef IMMERSED_BOUNDARY_H
#define IMMERSED_BOUNDARY_H

#include <bout/array.hxx>
#include <bout/bout_types.hxx>
#include <bout/field3d.hxx>
#include <bout/region.hxx>
#include <bout/utils.hxx>

/// Class for handling immersed boundary methods, meant for the
/// device wall.
class ImmersedBoundary {
public:
  ImmersedBoundary();
  enum class BoundCond {DIRICHLET, NEUMANN, SIZE};
  void SetBoundary(Field3D& f, const BoutReal bc, const BoundCond bc_type) const;
  BoutReal BoundaryNormalFlux(const Field3D& a, const Field3D& f, const Ind3D& i) const;
  bool IsInside(const Ind3D& ind) const;

private:
  /// Mask function that is 1 in the plasma, 0 in the wall
  Field3D bndry_mask;
  /// Mask function with ids into ghost cell data arrays.
  Field3D ghost_ids;
  /// Coordinate fields.
  Field3D R3;
  Field3D Z3;

  /// Ghost cell data arrays. #TODO: Add ghost_points in case needed. Useful for testing normals.
  int num_ghosts = 0;
  Matrix<BoutReal> image_inds;
  Matrix<BoutReal> ghost_points;
  Matrix<BoutReal> image_points;
  Matrix<BoutReal> bndry_points;
  Matrix<BoutReal> normals;
  Array<BoutReal> norm_dist;

  /// Image cell weights/ghost flag arrays.
  int num_weights = 0;
  Matrix<BoutReal> weights;
  Matrix<BoutReal> is_plasma;
  Array<int> all_plasma;

  //Boundary information.
  Field3D bound_ids;
  int num_bounds = 0;          // nb
  BoutReal s1 = 0;
  BoutReal s2 = 0;
  Array<BoutReal> bds;         // (nb)
  Array<BoutReal> bd_len;      // (nb)
  Matrix<BoutReal> mid_pts;    // (nb,2)
  Matrix<BoutReal> bnorms;     // (nb,2)
  Matrix<BoutReal> bbase_inds; // (nb,4): i0A,j0A,i0B,j0B
  Matrix<BoutReal> bweights;   // (nb, 2*nw): A then B

  /// TODO: Use these?
  /// Cell indices which are in the plasma.
  Region<Ind3D> plasma_region;
  /// Cell indices where a ghost point exists.
  Region<Ind3D> ghost_region;

  BoutReal GetGhostValue(const BoutReal image_val, const int gid,
                    const BoutReal bc, const BoundCond bc_type) const;
  BoutReal GetImageValue(Field3D& f, const int gid, const BoutReal bc_val,
                    const BoundCond bc_type) const;

  const std::string bc_exception = "Invalid boundary condition specified for immersed boundary.";

  // Solve 4x4 A x = b (A is copied by value; partial pivoting)
  // TODO: Clean up or find better solution...Thanks ChatGPT...
  template <typename T>
  static std::array<T,4> solve4x4(std::array<std::array<T,4>,4> A,
                                  const std::array<T,4> b) {
    int p[4] = {0,1,2,3};
    // LU with partial pivoting (Doolittle)
    for (int k = 0; k < 4; ++k) {
      // pivot
      int imax = k;
      T amax = std::abs(A[p[k]][k]);
      for (int i = k+1; i < 4; ++i) {
        T v = std::abs(A[p[i]][k]);
        if (v > amax) { amax = v; imax = i; }
      }
      std::swap(p[k], p[imax]);
      T pivot = A[p[k]][k];
      if (std::abs(pivot) == T(0)) { throw BoutException("Singular 4x4 Vandermonde"); }
      // eliminate
      for (int i = k+1; i < 4; ++i) {
        T m = A[p[i]][k] / pivot;
        A[p[i]][k] = m;
        for (int j = k+1; j < 4; ++j)
          A[p[i]][j] -= m * A[p[k]][j];
      }
    }
    // forward subst: Ly = b(p)
    T y[4];
    for (int i = 0; i < 4; ++i) {
      T sum = b[p[i]];
      for (int j = 0; j < i; ++j) sum -= A[p[i]][j] * y[j];
      y[i] = sum;
    }
    // back subst: Ux = y
    std::array<T,4> x;
    for (int i = 3; i >= 0; --i) {
      T sum = y[i];
      for (int j = i+1; j < 4; ++j) sum -= A[p[i]][j] * x[j];
      x[i] = sum / A[p[i]][i];
    }
    return x;
  }
};

//TODO Global for now.
extern ImmersedBoundary* immBdry;

#endif // IMMERSED_BOUNDARY_H
