
#include <bout/field.hxx>
#include <bout/globals.hxx>
#include <bout/constants.hxx>
#include <bout/mesh.hxx>
#include "../include/immersed_boundary.hxx"
using bout::globals::mesh;

//TODO: Make class if immersed == True? How to handle immersed flag.
//Can also run immersed BCs in each transform if fci.immersed == True.
ImmersedBoundary::ImmersedBoundary() {
  AUTO_TRACE();

  //TODO: Set up options based on // sheath bdry class?
  //Options& options = alloptions[name];
  xunit = sqrt(SI::qe*5.0/SI::Mp)/(SI::qe * 1.0 / SI::Mp);

  if (mesh->get(bndry_mask,     "in_mask") != 0 ||
      mesh->get(ghost_ids,     "ghost_id") != 0 ||
      mesh->get(num_ghosts,          "ng") != 0 ||
      mesh->get(num_weights,         "nw") != 0 ||
      mesh->get(image_inds,  "image_inds") != 0 ||
      mesh->get(image_points, "image_pts") != 0 ||
      mesh->get(bndry_points, "bndry_pts") != 0 ||
      mesh->get(normals,        "normals") != 0 ||
      mesh->get(norm_dist,    "norm_dist") != 0 || 
      mesh->get(is_plasma,    "is_plasma") != 0 ||
      mesh->get(weights,        "weights") != 0 ||
      mesh->get(R3,                   "R") != 0 ||
      mesh->get(Z3,                   "Z") != 0 ) {
        throw BoutException("Could not read immersed boundary values");
  }

  if (mesh->get(bound_ids,     "bound_id") != 0) {
    throw BoutException("Could not read bound_id values.");
  }
  if (mesh->get(num_bounds, "nb") != 0) {
    throw BoutException("Could not read nb values.");
  }
  if (mesh->get(mid_pts,     "mid_pts") != 0) {
    throw BoutException("Could not read mid_pts values.");
  }
  if (mesh->get(bnorms,     "bnorms") != 0) {
    throw BoutException("Could not read bnorms values.");
  }
  if (mesh->get(bd_len,     "bd_len") != 0) {
    throw BoutException("Could not read bd_len values.");
  }
  if (mesh->get(s1,     "s1") != 0) {
    throw BoutException("Could not read s1 values.");
  }
  if (mesh->get(s2,     "s2") != 0) {
    throw BoutException("Could not read s2 values.");
  }
  if (mesh->get(bweights,     "bweights") != 0) {
    throw BoutException("Could not read bweights values.");
  }
  if (mesh->get(bbase_inds, "base_inds") != 0) {
    throw BoutException("Could not read base_inds values.");
  }

  if (num_weights <= 0 or num_ghosts <= 0 or num_bounds <= 0) {
    throw BoutException("Invalid number of ghost cells or weights or bdy midpoints.");
  }

  //Set single flag if all neighbors are in plasma.
  all_plasma.reallocate(num_ghosts);
  for (size_t i = 0; i < all_plasma.size(); ++i) {
    bool is_all_plasma = true; //Assume true by default.
    for (size_t j = 0; j < num_weights; ++j) {
      is_all_plasma &= static_cast<bool>(is_plasma(i,j)); //TODO how to access at i earlier?
    }
    all_plasma[i] = is_all_plasma;
  }
}

bool ImmersedBoundary::IsInside(const Ind3D& ind) const {
  return static_cast<bool>(bndry_mask(ind.x(),0,ind.z()));
}

BoutReal ImmersedBoundary::BoundaryNormalFlux(const Field3D& a, const Field3D& f, const Ind3D& i) const {
  const auto x = i.x();
  const auto y = i.y();
  const auto z = i.z();

  const int b = static_cast<int>(std::lround(bound_ids(x,0,z)));
  if (b < 0) {
    return 0.0;
  }

  BoutReal dfdn = 0.0;
  const auto bc_type = BoundCond::NEUMANN;
  if (bc_type == BoundCond::NEUMANN) {
    dfdn = 2.0; // df/dn at boundary, use negative for outward normal here. TODO: Get from input with correct sign?.
    const auto bflux = -a[i] * dfdn * bd_len[b] * 1.0;
    return bflux; //TODO: Dont need to duplicate this logic?
  }

  // --- Base indices for A and B (lower-left cell center of bilinear stencil)
  //TODO: Get ints reading in correctly...
  const int ai0 = static_cast<int>(bbase_inds(b, 0));
  const int aj0 = static_cast<int>(bbase_inds(b, 1));
  const int bi0 = static_cast<int>(bbase_inds(b, 2));
  const int bj0 = static_cast<int>(bbase_inds(b, 3));

  // Weights (your ordering: w00,w01,w10,w11)
  const BoutReal wA00 = bweights(b, 0);
  const BoutReal wA01 = bweights(b, 1);
  const BoutReal wA10 = bweights(b, 2);
  const BoutReal wA11 = bweights(b, 3);

  const BoutReal wB00 = bweights(b, 4);
  const BoutReal wB01 = bweights(b, 5);
  const BoutReal wB10 = bweights(b, 6);
  const BoutReal wB11 = bweights(b, 7);

  // Bilinear interpolation using the correct local stencil
  const BoutReal fA = wA00*f(ai0,1,aj0) + wA01*f(ai0,1,aj0+1) + wA10*f(ai0+1,1,aj0) + wA11*f(ai0+1,1,aj0+1);
  const BoutReal fB = wB00*f(bi0,1,bj0) + wB01*f(bi0,1,bj0+1) + wB10*f(bi0+1,1,bj0) + wB11*f(bi0+1,1,bj0+1);

  // Boundary value and one-sided derivative from J&C.
  const BoutReal fb = 1.0/3.0; //.09; //TODO: Sample at midpoint for J&C.

  const BoutReal term1 = (s2/s1) * (fb - fA);
  const BoutReal term2 = (s1/s2) * (fb - fB);

  dfdn = (term1 - term2) / (s2 - s1);

  //Flux value.
  //Use negative because normals oriented inward. Positive default for outward normals.
  const auto bflux = -a[i] * bd_len[b] * 1.0 * dfdn;
  return bflux; //TODO: Evaluate a at the boundary (const now so OK). (And move out of IB to div_ops!?)
}

// Calculate image value from nearby grid points. Note weights are defined as
// w00, w01, w10, w11 with indices (x,z). All other values follow from there.
BoutReal ImmersedBoundary::GetImageValue(Field3D& f, const int gid,
            const BoutReal bc_val, const BoundCond bc_type) const {
  // Get nearby vals to image from floating point index.
  int indx = static_cast<int>(image_inds(gid,0));
  int indz = static_cast<int>(image_inds(gid,1));
  //TODO: Use BOUT style arrays?
  auto node_vals = std::array<BoutReal, 4>{f(indx,1,indz), f(indx,1,indz+1), //TODO: Need indy for ghost cells in y? Same below.
                                           f(indx+1,1,indz), f(indx+1,1,indz+1)};

  BoutReal image_val = 0.0;
  // If all nearby nodes in plasma just add weights.
  if (all_plasma[gid]) {
    //TODO: Get weights[gid] first? Same with plasma_flags below...
    //TODO: Can also combine code and solve with Vandermonde just fine. Set up beforehand to do a dot product.
    for (size_t i = 0; i < num_weights; ++i) {
      image_val += weights(gid,i)*node_vals[i];
    }
  }
  // If some nearby points are ghost cells.
  else {
    //TODO: Deal with not knowing num_weights?
    std::array<std::array<BoutReal, 4>, 4> vandMat{};
    // Get R,Z values of nearby cells. //TODO: Use Ind3D and xp(), etc...?
    auto nodes_x = std::array<BoutReal, 4>{R3(indx,0,indz),   R3(indx,0,indz+1),
                                          R3(indx+1,0,indz), R3(indx+1,0,indz+1)};
    auto nodes_z = std::array<BoutReal, 4>{Z3(indx,0,indz),   Z3(indx,0,indz+1),
                                          Z3(indx+1,0,indz), Z3(indx+1,0,indz+1)};
    // Get indices to nearby ghost cell data.
    auto node_gids = std::array<BoutReal, 4>{ghost_ids(indx,0,indz), ghost_ids(indx,0,indz+1),
                                             ghost_ids(indx+1,0,indz), ghost_ids(indx+1,0,indz+1)};

    for (size_t i = 0; i < num_weights; ++i) {
      int igid = static_cast<int>(std::lround(node_gids[i])); //TODO Template to allow for reading ints and remove all lrounds...
      if (gid == igid) { //If primary ghost cell.
        auto node_gid = igid;
        auto xB = bndry_points(node_gid,0);
        auto zB = bndry_points(node_gid,1);
        auto xN = normals(node_gid,0);
        auto zN = normals(node_gid,1);
        node_vals[i] = bc_val; // Just change the node val directly.
        switch (bc_type) {
          case BoundCond::DIRICHLET:
            vandMat[i] = std::array<BoutReal, 4>{xB*zB, xB, zB, 1.0};
            break;
          case BoundCond::NEUMANN:
            node_vals[i] = -node_vals[i]; //TODO: Fix neg sign? Same as below...
            vandMat[i] = std::array<BoutReal, 4>{xB*zN + zB*xN, xN, zN, 0};
              break;
          default:
            throw BoutException(bc_exception);
        }
      } else {
        auto x = nodes_x[i];
        auto z = nodes_z[i];
        vandMat[i] = std::array<BoutReal, 4>{x*z, x, z, 1.0};
      }
    }

    // Perform 4x4 matrix solve.
    // TODO: Setup everything once for each ghost point so quick to solve. And use PETSC?
    // ChatGPT seems to think it can be set up so only a dot product is needed each timestep.
    auto c = solve4x4(vandMat, node_vals);

    auto xI = image_points(gid,0);
    auto zI = image_points(gid,1);
    // Set image value from solved coefficients.
    //TODO: Double check things seem ok for both BCs.
    image_val = c[0]*(xI*zI) + c[1]*xI + c[2]*zI + c[3];
  }

  return image_val;
}

BoutReal ImmersedBoundary::GetGhostValue(const BoutReal image_val, const int gid,
                          const BoutReal bc_val, const BoundCond bc_type) const {
  switch (bc_type) {
    case BoundCond::DIRICHLET:
      return 2*bc_val - image_val;
    case BoundCond::NEUMANN:
      return image_val - 2.0*norm_dist[gid]*(-bc_val); //TODO: Fix neg sign?
    default:
        throw BoutException(bc_exception);
  }
}

void ImmersedBoundary::SetBoundary(Field3D& f, const BoutReal bc_val,
                                   const BoundCond bc_type) const {
  // Do a few Gauss–Seidel-style sweeps so ghost cells that depend on other
  // ghost cells (via interpolation) can converge. //TODO:Get mpi working.
  constexpr int max_gs_iters = 12;

  for (int it = 0; it < max_gs_iters; ++it) { //TODO: Is this equivalent to GS iteration?
    // NOTE: This is "Gauss–Seidel-like" because we update f in-place as we sweep.
    BOUT_FOR_SERIAL(i, f.getRegion("RGN_NOBNDRY")) {
      const auto x = i.x();
      const auto y = i.y();
      const auto z = i.z();

      const auto gid = ghost_ids(x, 0, z);
      if (gid >= 0) {
        if (it == 0 || //TODO: Clean up this if check?
           (it > 0 && (is_plasma(gid,0) + is_plasma(gid,1) + is_plasma(gid,2) + is_plasma(gid,3) < 3))) {
            const auto image_val = GetImageValue(f, gid, bc_val, bc_type);
            const auto ghost_val = GetGhostValue(image_val, gid, bc_val, bc_type);
            f(x, y, z) = ghost_val;
        }
      } else if (bndry_mask(x, 0, z) == 0) {
          f(x, y, z) = 0.0;
      }
    }
  }
}