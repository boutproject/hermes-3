#include <../include/immersed_boundary.hxx>

#include <regex>

#include <bout/field.hxx>
#include <bout/globals.hxx>
#include <bout/constants.hxx>
#include <bout/mesh.hxx>
using bout::globals::mesh;

ImmersedBoundary::ImmersedBoundary() {
  AUTO_TRACE();

  ASSERT0(mesh->get(bndry_mask,     "in_mask") == 0);
  ASSERT0(mesh->get(ghost_ids,     "ghost_id") == 0);
  ASSERT0(mesh->get(num_ghosts,          "ng") == 0);
  ASSERT0(mesh->get(num_weights,         "nw") == 0);
  ASSERT0(mesh->get(image_inds,  "image_inds") == 0);
  ASSERT0(mesh->get(image_points, "image_pts") == 0);
  ASSERT0(mesh->get(bndry_points, "bndry_pts") == 0);
  ASSERT0(mesh->get(normals,        "normals") == 0);
  ASSERT0(mesh->get(norm_dist,    "norm_dist") == 0);
  ASSERT0(mesh->get(is_plasma,    "is_plasma") == 0);
  ASSERT0(mesh->get(weights,        "weights") == 0);
  ASSERT0(mesh->get(R3,                   "R") == 0);
  ASSERT0(mesh->get(Z3,                   "Z") == 0);
  ASSERT0(mesh->get(bound_ids,     "bound_id") == 0);
  ASSERT0(mesh->get(num_bounds,          "nb") == 0);
  ASSERT0(mesh->get(mid_pts,        "mid_pts") == 0);
  ASSERT0(mesh->get(bnorms,          "bnorms") == 0);
  ASSERT0(mesh->get(bd_len,          "bd_len") == 0);
  ASSERT0(mesh->get(s1,                  "s1") == 0);
  ASSERT0(mesh->get(s2,                  "s2") == 0);
  ASSERT0(mesh->get(bweights,      "bweights") == 0);
  ASSERT0(mesh->get(bbase_inds,   "base_inds") == 0);

  if (num_weights <= 0 or num_ghosts <= 0 or num_bounds <= 0) {
    throw BoutException("Invalid number of ghost cells or weights or cut cells.");
  }

  //Set useful flags for info about image nodes being plasma cells or not.
  all_plasma.reallocate(num_ghosts);
  ghost_count.reallocate(num_ghosts);
  for (size_t i = 0; i < num_ghosts; ++i) {
    bool is_all_plasma = true; //Assume true by default.
    ghost_count[i] = num_weights;
    for (size_t j = 0; j < num_weights; ++j) {
      const auto is_plasma_node = static_cast<bool>(is_plasma(i,j));
      is_all_plasma &= is_plasma_node;
      ghost_count[i] -= is_plasma_node;
    }
    all_plasma[i] = is_all_plasma;
  }

  //Set up new boundary and plasma regions.
  Region<Ind3D>::RegionIndices xbndry_indices;
  Region<Ind3D>::RegionIndices zbndry_indices;
  Region<Ind3D>::RegionIndices nobndry_indices;
  Region<Ind3D>::RegionIndices nobndrycut_indices;
  Region<Ind3D>::RegionIndices bndry_indices;
  Region<Ind3D>::RegionIndices gst_indices;
  BOUT_FOR(i, bndry_mask.getRegion("RGN_NOBNDRY")) {
    //Ignore guard cells in y for which cells to evolve.
    if (i.y() < mesh->ystart || i.y() > mesh->yend) {continue;}
    if (IsInside(i)) {
      nobndry_indices.push_back(i);
      if (IsCutCell(i)) {nobndrycut_indices.push_back(i);}
    }
    else {
      if (IsGhost(i)) {
        gst_indices.push_back(i);
        //Note, dont if/else here because need both indices in both regions sometimes.
        if (IsInside(i.xp())) {xbndry_indices.push_back(i);}
        if (IsInside(i.zp())) {zbndry_indices.push_back(i);}
      }
      bndry_indices.push_back(i);
    }
  }

  //IMM_BNDRY_TODO: Need upsert just to replace RGN_NOBNDRY if considered worthwhile.
  //FV operators work on +x/+z faces only so subset of ghosts here.
  mesh->upsertRegion3D("RGN_dagp_fv_xbndry",   Region<Ind3D>(xbndry_indices));
  mesh->upsertRegion3D("RGN_dagp_fv_zbndry",   Region<Ind3D>(zbndry_indices));
  //Ghost indices.
  mesh->upsertRegion3D("RGN_IMM_BNDRY_GST",    Region<Ind3D>(gst_indices));
  //Plasma indices.
  mesh->upsertRegion3D("RGN_NO_IMM_BNDRY",     Region<Ind3D>(nobndry_indices));
  mesh->upsertRegion3D("RGN_NO_IMM_BNDRY_CUT", Region<Ind3D>(nobndrycut_indices));
  //Non-plasma indices (including ghosts).
  mesh->upsertRegion3D("RGN_IMM_BNDRY",        Region<Ind3D>(bndry_indices));
}

bool ImmersedBoundary::IsGhost(const Ind3D& ind) const {
  return ghost_ids(ind.x(),ind.y(),ind.z()) >= 0;
}

bool ImmersedBoundary::IsInside(const Ind3D& ind) const {
  return static_cast<bool>(bndry_mask(ind.x(),ind.y(),ind.z()));
}

bool ImmersedBoundary::IsCutCell(const Ind3D& ind) const {
  return bound_ids(ind.x(),ind.y(),ind.z()) >= 0;
}

ImmersedBoundary::BC_Info ImmersedBoundary::ReadBC(const std::string& bc_info) const {
  //IMM_BNDRY_TODO, the way this is done should probably pull from the boundary factory code and be shared.
  //Match input BC expression. Allows for neumann/dirichlet, scientific notation and various brackets.
  static const std::regex bc_value_regex(
    R"(\s*(neumann|dirichlet)\s*[\(\[\{\<]\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)\s*[\)\]\}\>]\s*)"
  );
  static const std::unordered_map<std::string, BoundCond> str_to_bc {
    { "neumann",  BoundCond::NEUMANN },
    { "dirichlet", BoundCond::DIRICHLET }
  };

  std::smatch match;
  if (std::regex_match(bc_info, match, bc_value_regex)) {
    const std::string bc_type = match[1];
    const BoutReal value = std::stod(match[2]);

    auto it = str_to_bc.find(bc_type);
    if (it == str_to_bc.end()) {
      throw BoutException("Unknown BC type for immersed boundary: '" + bc_type + "'");
    }

    const auto& bc_info = BC_Info(it->second, value);
    return bc_info;
  } else {
    throw BoutException("Unknown BC input for immersed boundary: '" + bc_info + "'");
  }
}

void ImmersedBoundary::FieldSetup(Field3D& f) {
  const auto& bc_type = Options::root()[f.name][bc_key]
        .doc("Boundary condition to use at immersed boundary wall.")
        .withDefault("neumann(0.0)"); //Default to no flux conditions.

  f.setRegion("RGN_NO_IMM_BNDRY"); //Default region to plasma cells.

  //Load BC info and clear out non-plasma cells.
  const auto& bc_info = ReadBC(bc_type);
  bc_map[f.name] = bc_info;
  BOUT_FOR(i, f.getRegion("RGN_IMM_BNDRY")) {
    f[i] = 0.0;
    ddt(f)[i] = 0.0;
  }
}

//IMM_BNDRY_TODO: Interpolate dfdn from dfdn at segment endpoints (bdy/cell intersections).
//IMM_BNDRY_TODO: Also need fb and a at boundary approximation.
void ImmersedBoundary::ComputeBoundaryFluxes(const Field3D& a, const Field3D& f, Field3D& result) const {
  //Get bc info from field first.
  const auto bc_type = bc_map.at(f.name).first;
  const auto bc_val  = bc_map.at(f.name).second;

  //Loop over cut-cells w/ centers in plasma and calculate flux through wall boundary.
  BOUT_FOR(i, f.getRegion("RGN_NO_IMM_BNDRY_CUT")) {
    const auto b = bound_ids[i];

    BoutReal dfdn = 0.0;
    if (bc_type == BoundCond::NEUMANN) {
      dfdn = bc_val; //Note, using outward normal derivative here.
    } else if (bc_type == BoundCond::DIRICHLET) {
      //Perform bilinear interpolation for points along the boundary normal, and 
      //calculate normal derivative at boundary.
      const int ai0 = static_cast<int>(bbase_inds(b, 0));
      const int aj0 = static_cast<int>(bbase_inds(b, 1));
      const int bi0 = static_cast<int>(bbase_inds(b, 2));
      const int bj0 = static_cast<int>(bbase_inds(b, 3));

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

      // Boundary value and one-sided derivative from J&C '98.
      const BoutReal fb = bc_val;

      const BoutReal term1 = (s2/s1) * (fb - fA);
      const BoutReal term2 = (s1/s2) * (fb - fB);

      dfdn = (term1 - term2) / (s2 - s1);
    } else {
       throw BoutException(bc_exception);
    }

    //Flux assumes normal is outward for dfdn.
    const BoutReal depthFac = 1.0; //IMM_BNDRY_TODO: Use toroidal curvature (at cell center?).
    const auto bdy_flux = -a[i] * bd_len[b] * depthFac * dfdn;
    result[i] -= bdy_flux;
  }
}

// Calculate image value from nearby grid points. Note weights are defined as
// w00, w01, w10, w11 with indices (x,z). All other values follow from there.
BoutReal ImmersedBoundary::GetImageValue(Field3D& f, const int gid,
            const BoutReal bc_val, const BoundCond bc_type) const {
  // Get nearby vals to image from floating point index.
  int indx = static_cast<int>(image_inds(gid,0));
  int indz = static_cast<int>(image_inds(gid,1));
  //TODO: Need indy for ghost cells in y? Same below. Or do everything in 2d?
  //TODO: Use num_weights instead of 4, same for below where 4s used.
  auto node_vals = std::array<BoutReal, 4>{f(indx,1,indz), f(indx,1,indz+1), 
                                           f(indx+1,1,indz), f(indx+1,1,indz+1)};

  BoutReal image_val = 0.0;
  // If all nearby nodes in plasma just add weights.
  if (all_plasma[gid]) {
    for (size_t i = 0; i < num_weights; ++i) {
      image_val += weights(gid,i)*node_vals[i];
    }
  }
  // If some nearby points are ghost cells.
  else {
    std::array<std::array<BoutReal, 4>, 4> vandMat{};
    // Get R,Z values of nearby cells.
    auto nodes_x = std::array<BoutReal, 4>{R3(indx,1,indz),   R3(indx,1,indz+1),
                                           R3(indx+1,1,indz), R3(indx+1,1,indz+1)};
    auto nodes_z = std::array<BoutReal, 4>{Z3(indx,1,indz),   Z3(indx,1,indz+1),
                                           Z3(indx+1,1,indz), Z3(indx+1,1,indz+1)};
    // Get indices to nearby ghost cell data.
    auto node_gids = std::array<int, 4>{static_cast<int>(std::lround(ghost_ids(indx,1,indz))),
                                        static_cast<int>(std::lround(ghost_ids(indx,1,indz+1))),
                                        static_cast<int>(std::lround(ghost_ids(indx+1,1,indz))),
                                        static_cast<int>(std::lround(ghost_ids(indx+1,1,indz+1)))};

    for (size_t i = 0; i < num_weights; ++i) {
      int igid = node_gids[i];

      //If primary ghost cell || singular secondary ghost cell.
      if (gid == igid || (igid >= 0 && ghost_count[igid] == 1)) {
        auto node_gid = igid;
        auto xB = bndry_points(node_gid,0);
        auto zB = bndry_points(node_gid,1);
        auto xN = normals(node_gid,0);
        auto zN = normals(node_gid,1);
        switch (bc_type) {
          case BoundCond::DIRICHLET:
            node_vals[i] = bc_val;
            vandMat[i] = std::array<BoutReal, 4>{xB*zB, xB, zB, 1.0};
            break;
          case BoundCond::NEUMANN:
            node_vals[i] = -bc_val;
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

    //Perform 4x4 matrix solve.
    //TODO: Combine code to always use Vandermonde, and solve with Vandermonde beforehand to do a dot product.
    auto c = solve4x4(vandMat, node_vals);

    auto xI = image_points(gid,0);
    auto zI = image_points(gid,1);

    //Set image value from solved coefficients.
    image_val = c[0]*(xI*zI) + c[1]*xI + c[2]*zI + c[3];
  }

  return image_val;
}

BoutReal ImmersedBoundary::GetGhostValue(const BoutReal image_val, const int gid,
                          const BoutReal bc_val, const BoundCond bc_type) const {
  switch (bc_type) {
    case BoundCond::DIRICHLET:
      return 2.0*bc_val - image_val;
    case BoundCond::NEUMANN:
      return image_val - 2.0*norm_dist[gid]*(-bc_val);
    default:
        throw BoutException(bc_exception);
  }
}

//IMM_BNDRY_TODO: Get mpi working.
void ImmersedBoundary::SetBoundary(Field3D& f) {
  const auto bc_type = bc_map.at(f.name).first;
  const auto bc_val  = bc_map.at(f.name).second;

  //Multiple ghost node case requires Gauss-Seidel convergence due to coupled initial guesses.
  constexpr int max_gs_iters = 12; //TODO: Make a vector to store old values, reset to 0 at start, and compare to stop sweeping.
  for (int it = 0; it < max_gs_iters; ++it) {
    //NOTE: This is "Gauss–Seidel-like" because we update f in-place as we sweep.
    BOUT_FOR(i, f.getRegion("RGN_IMM_BNDRY_GST")) {
      const auto gid = ghost_ids[i];
      if (gid >= 0) {
        //If first iteration or more than 1 ghost (GS iteration to converge).
        if (it == 0 || (it > 0 && ghost_count[gid] > 1)) {
          const auto image_val = GetImageValue(f, gid, bc_val, bc_type);
          const auto ghost_val = GetGhostValue(image_val, gid, bc_val, bc_type);
          f[i] = ghost_val;
        }
      }
    }
  }
}