
#include <bout/field.hxx>
#include <bout/globals.hxx>
#include <bout/mesh.hxx>
#include "../include/immersed_boundary.hxx"
using bout::globals::mesh;

//TODO: Make class if immersed == True? How to handle immersed flag.
//Can also run immersed BCs in each transform if fci.immersed == True.
ImmersedBoundary::ImmersedBoundary() {
  AUTO_TRACE();

  //TODO: Set up options based on // sheath bdry class?
  //Options& options = alloptions[name];

  if (mesh->get(bndry_mask,     "in_mask") != 0 ||
      mesh->get(ghost_ids,     "ghost_id") != 0 ||
      mesh->get(num_ghosts,          "ng") != 0 ||
      mesh->get(num_weights,         "nw") != 0 ||
      mesh->get(image_inds,  "image_inds") != 0 ||
      //mesh->get(ghost_points, "ghost_pts") != 0 ||
      mesh->get(image_points, "image_pts") != 0 ||
      mesh->get(bndry_points, "bndry_pts") != 0 ||
      mesh->get(normals,        "normals") != 0 ||
      mesh->get(norm_dist,    "norm_dist") != 0 || 
      mesh->get(is_plasma,    "is_plasma") != 0 ||
      mesh->get(weights,        "weights") != 0 ||
      mesh->get(R3,                   "R") != 0 ||
      mesh->get(Z3,                   "Z") != 0) {
        throw BoutException("Could not read immersed boundary values");
  }

  if (num_weights <= 0 or num_ghosts <= 0) {
    throw BoutException("Invalid number of ghost cells or weights.");
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
  return static_cast<bool>(bndry_mask(ind.x(), 0, ind.z()));
}

void ImmersedBoundary::ResetGhosts(Field3D& f) const {
  floor(f, 0.0);
}

// Calculate image value from nearby grid points. Note weights are defined as
// w00, w01, w10, w11 with indices (x,z). All other values follow from there.
BoutReal ImmersedBoundary::GetImageValue(Field3D& f, const int gid,
            const BoutReal bc_val, const BoundCond bc_type) const {
  // Get nearby vals to image from floating point index.
  int indx = static_cast<int>(image_inds(gid,0));
  int indz = static_cast<int>(image_inds(gid,1));
  //TODO: Use BOUT style arrays?
  auto node_vals = std::array<BoutReal, 4>{f(indx,0,indz), f(indx,0,indz+1), //TODO: Need indy for ghost cells in y? Same below.
                                           f(indx+1,0,indz), f(indx+1,0,indz+1)};

  BoutReal image_val = 0.0;
  // If all nearby nodes in plasma just add weights.
  if (all_plasma[gid]) {
    //TODO: Get weights[gid] first? Same with plasma_flags below...
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
      if (!is_plasma(gid,i)) {
        auto node_gid = node_gids[i];
        auto xB = bndry_points(node_gid,0);
        auto zB = bndry_points(node_gid,1);
        auto xN = normals(node_gid,0);
        auto zN = normals(node_gid,1);
        node_vals[i] = bc_val; // Just change the node val directly.
        switch (bc_type) {
          case BoundCond::DIRICHLET:
            vandMat[i] = std::array<BoutReal, 4>{xB*zB, zB, xB, 1};
            break;
          case BoundCond::NEUMANN:
            vandMat[i] = std::array<BoutReal, 4>{xN*zB + zN*xB, zN, xN, 0};
              break;
          default:
            throw BoutException(bc_exception);
        }
      }
      else {
        auto x = nodes_x[i];
        auto z = nodes_z[i];
        vandMat[i] = std::array<BoutReal, 4>{x*z, z, x, 1};
      }
    }

    // Perform 2x2 matrix solve.
    // TODO: Setup everything once for each ghost point so quick to solve. And use PETSC?
    // ChatGPT seems to think it can be set up so only a dot product is needed each timestep.
    auto c = solve4x4(vandMat, node_vals);

    auto xI = image_points(gid,0);
    auto zI = image_points(gid,1);
    // Set image value from solved coefficients.
    //TODO: Double check things seem ok for both BCs.
    image_val = c[0]*(xI*zI) + c[1]*zI + c[2]*xI + c[3];
  }

  return image_val;
}

BoutReal ImmersedBoundary::GetGhostValue(const BoutReal image_val, const int gid,
                          const BoutReal bc_val, const BoundCond bc_type) const {
  switch (bc_type) {
    case BoundCond::DIRICHLET:
      return 2*bc_val - image_val;
    case BoundCond::NEUMANN:
      return image_val - 2*norm_dist[gid]*bc_val;
    default:
        throw BoutException(bc_exception);
  }
}

void ImmersedBoundary::SetBoundary(Field3D& f, const BoutReal bc_val,
                                     const BoundCond bc_type) const {
  //TODO: Need to change if 3D...Can ignore y ghosts or not?
  BOUT_FOR_SERIAL(i, f.getRegion("RGN_NOBNDRY")) {
    auto x = i.x();
    auto y = i.y();
    auto z = i.z();
    auto gid = ghost_ids(x,0,z);
    if (gid >= 0) {
        //std::cout << "gid: " << gid << std::endl;
        auto image_val = GetImageValue(f, gid, bc_val, bc_type);
        auto ghost_val = GetGhostValue(image_val, gid, bc_val, bc_type);
        f(x,y,z) = ghost_val;
    } else {
      if (bndry_mask(x,0,z) == 0) {
        f(x,y,z) = 0.0;
      }
    }
  }
}