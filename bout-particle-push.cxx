#include "bout/bout.hxx"
#include "bout/field2d.hxx"
#include "bout/output.hxx"
#include "bout/petsclib.hxx"
#include <bout/field_factory.hxx>
#include <fmt/core.h>
#include <iostream>
#include <neso_particles.hpp>
#include <neso_particles/external_interfaces/petsc/petsc_interface.hpp>
#include <neso_particles/typedefs.hpp>
#include <netcdf>
#include <petscsystypes.h>
#include <petscviewerhdf5.h>
#include <string>
#include <vector>

#ifndef NESO_PARTICLES_PETSC
static_assert(false, "NESO-Particles was installed without PETSc support.");
#else

template <typename T, typename U>
inline void ASSERT_EQ(T t, U u) {
  NESOASSERT(t == u, "A check failed.");
}

void collect_unique_points(std::vector<double>& global_Z_vertices_buffer,
                           std::vector<double>& global_R_vertices_buffer, int& N_unique,
                           const double& zero,
                           std::vector<double>& global_Z_hypnotoad_vertices,
                           std::vector<double>& global_R_hypnotoad_vertices) {
  bool unique;
  int N_nonunique_vertices = global_Z_hypnotoad_vertices.size();
  for (int iv = 0; iv < N_nonunique_vertices; iv++) {
    // assume the point
    // iv = (global_R_hypnotoad_vertices.at(iv),global_Z_hypnotoad_vertices.at(iv))
    // is unique
    unique = true;
    // check if the point is unique, by comparing the the existing N_unique points
    for (int iunique = 0; iunique < N_unique; iunique++) {
      if (abs(global_Z_hypnotoad_vertices.at(iv) - global_Z_vertices_buffer.at(iunique))
              < zero
          && abs(global_R_hypnotoad_vertices.at(iv)
                 - global_R_vertices_buffer.at(iunique))
                 < zero) {
        unique = false;
        // we have determined that the point is not unique
      }
    }
    if (unique) {
      // add the point and increment N_unique
      global_Z_vertices_buffer.at(N_unique) = global_Z_hypnotoad_vertices.at(iv);
      global_R_vertices_buffer.at(N_unique) = global_R_hypnotoad_vertices.at(iv);
      N_unique++;
    }
  }
}

void RZ_to_ivertex_vector(Field2D& ivertex_corners,
                          std::vector<double>& global_Z_vertices,
                          std::vector<double>& global_R_vertices, const double& zero,
                          Mesh*& bout_mesh, Field2D& Rxy_corners, Field2D& Zxy_corners) {
  int Nvertex = global_Z_vertices.size();
  bool index_found;
  // loop over points that are not guard cells
  for (int ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (int iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      index_found = false;
      for (int iv = 0; iv < Nvertex; iv++) {
        if (abs(global_Z_vertices.at(iv) - Zxy_corners(ix, iy)) < zero
            && abs(global_R_vertices.at(iv) - Rxy_corners(ix, iy)) < zero) {
          ivertex_corners(ix, iy) = iv;
          index_found = true;
          // we have matched an iv global index to a (R,Z) from Hypnotoad
        }
      }
      // exit if we fail to find a match
      NESOASSERT(index_found,
                 fmt::format("ivertex not found for ix = {} iy = {}", ix, iy));
    }
  }
}

void load_vertex_information_from_netcdf(int& Nvertex,
                                         std::vector<double>& global_vertex_R,
                                         std::vector<double>& global_vertex_Z) {
  // read data from netcdf for global vertices in mesh
  // Open the NetCDF file in read-only mode
  const std::string filename = Options::root()["mesh"]["file"];
  netCDF::NcFile dataFile(filename, netCDF::NcFile::read);

  // Get the variable
  std::string varName = "global_vertex_list_R";
  netCDF::NcVar dataVar = dataFile.getVar(varName);
  NESOASSERT(!dataVar.isNull(), fmt::format("Variable {} not found in file.", varName));
  std::vector<netCDF::NcDim> dims = dataVar.getDims();
  size_t nvertices = dims[0].getSize();

  // Read the data into a vector
  std::vector<double> global_vertex_list_R(nvertices);
  dataVar.getVar(global_vertex_list_R.data());

  // Get the variable
  varName = "global_vertex_list_Z";
  dataVar = dataFile.getVar(varName);
  NESOASSERT(!dataVar.isNull(), fmt::format("Variable {} not found in file.", varName));
  // Read the data into a vector
  std::vector<double> global_vertex_list_Z(nvertices);
  dataVar.getVar(global_vertex_list_Z.data());
  dataFile.close();
  // assign data to output variables
  Nvertex = nvertices;
  global_vertex_R = global_vertex_list_R;
  global_vertex_Z = global_vertex_list_Z;
}

std::vector<PetscInt> cells_definition_from_RZ_ivertex(
    std::vector<PetscInt>& cells, Mesh*& bout_mesh, Field2D& Rxy_lower_left_corners,
    Field2D& Rxy_lower_right_corners, Field2D& Rxy_upper_right_corners,
    Field2D& Rxy_upper_left_corners, Field2D& Zxy_lower_left_corners,
    Field2D& Zxy_lower_right_corners, Field2D& Zxy_upper_right_corners,
    Field2D& Zxy_upper_left_corners, Field2D& ivertex_lower_left_corners,
    Field2D& ivertex_lower_right_corners, Field2D& ivertex_upper_right_corners,
    Field2D& ivertex_upper_left_corners) {
  std::vector<PetscReal> Z_vertices(4);
  std::vector<PetscReal> R_vertices(4);
  std::vector<PetscReal> theta_vertices(4);
  std::vector<PetscInt> i_vertices(4);
  std::vector<PetscInt> sort_indices(4);
  PetscReal ZZ;
  PetscReal ZZmid;
  PetscReal RR;
  PetscReal RRmid;
  // local number of x cells, excluding guards
  int Nx = bout_mesh->xend - bout_mesh->xstart + 1;
  // local number of y cells, excluding guards
  int Ny = bout_mesh->yend - bout_mesh->ystart + 1;
  PetscInt num_cells_owned = Nx * Ny;
  // std::vector<PetscInt> cells;
  cells.reserve(num_cells_owned * 4);
  // We are careful to list the vertices in counter clock-wise order.
  for (PetscInt ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (PetscInt iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      // collect data from Hypnotoad arrays
      i_vertices[0] = static_cast<int>(std::lround(ivertex_lower_left_corners(ix, iy)));
      i_vertices[1] = static_cast<int>(std::lround(ivertex_lower_right_corners(ix, iy)));
      i_vertices[2] = static_cast<int>(std::lround(ivertex_upper_right_corners(ix, iy)));
      i_vertices[3] = static_cast<int>(std::lround(ivertex_upper_left_corners(ix, iy)));

      R_vertices[0] = Rxy_lower_left_corners(ix, iy);
      R_vertices[1] = Rxy_lower_right_corners(ix, iy);
      R_vertices[2] = Rxy_upper_right_corners(ix, iy);
      R_vertices[3] = Rxy_upper_left_corners(ix, iy);

      Z_vertices[0] = Zxy_lower_left_corners(ix, iy);
      Z_vertices[1] = Zxy_lower_right_corners(ix, iy);
      Z_vertices[2] = Zxy_upper_right_corners(ix, iy);
      Z_vertices[3] = Zxy_upper_left_corners(ix, iy);
      // get the mean values of R, Z
      RRmid = 0.0;
      ZZmid = 0.0;
      for (PetscInt iv = 0; iv < 4; iv++) {
        RRmid += R_vertices[iv];
        ZZmid += Z_vertices[iv];
      }
      RRmid /= 4.0;
      ZZmid /= 4.0;
      // get the angle subtended from the centre of the cell to each vertex
      for (PetscInt iv = 0; iv < 4; iv++) {
        RR = R_vertices[iv] - RRmid;
        ZZ = Z_vertices[iv] - ZZmid;
        theta_vertices[iv] = std::atan2(ZZ, RR);
      }
      // get the indices that sort the vertices in ascending order of theta
      for (PetscInt iv = 0; iv < 4; ++iv) {
        sort_indices[iv] = iv;
      }
      std::sort(sort_indices.begin(), sort_indices.end(),
                [&theta_vertices](int i, int j) {
                  return theta_vertices[i] < theta_vertices[j];
                });
      // fill cells using the sorted indices
      for (PetscInt iv = 0; iv < 4; ++iv) {
        cells.push_back(i_vertices[sort_indices[iv]]);
      }
    }
  }
  return cells;
}

using namespace NESO::Particles;

int main(int argc, char** argv) {
  // initialise_mpi(&argc, &argv);
  // attempt to call BOUT to
  // get information from a BOUT
  // mesh object
  // N.B. Comment the next three lines
  // to permit compilation as is
  BoutInitialise(argc, argv);
  Mesh* bout_mesh = Mesh::create(&Options::root()["mesh"]);
  bool use_cxx_ivertex = Options::root()["mesh"]["use_cxx_ivertex"].withDefault(false);
  std::string dmplex_name =
      Options::root()["mesh"]["dmplex_name"].withDefault("hypnotoad_dmplex_mesh");
  std::string dmplex_h5_filename =
      Options::root()["mesh"]["dmplex_h5_filename"].withDefault(
          "hypnotoad_dmplex_mesh_output.h5");
  output << fmt::format("Using option use_cxx_ivertex = {}", use_cxx_ivertex)
         << std::endl;
  bout_mesh->load();
  Coordinates* coord = bout_mesh->getCoordinates();
  Field2D Rxy_lower_left_corners;
  Field2D Rxy_lower_right_corners;
  Field2D Rxy_upper_right_corners;
  Field2D Rxy_upper_left_corners;
  Field2D Zxy_lower_left_corners;
  Field2D Zxy_lower_right_corners;
  Field2D Zxy_upper_right_corners;
  Field2D Zxy_upper_left_corners;
  // mesh->get(ivertex, "ivertex_lower_left_corners");
  bout_mesh->get(Rxy_lower_left_corners, "Rxy_corners");
  bout_mesh->get(Rxy_lower_right_corners, "Rxy_lower_right_corners");
  bout_mesh->get(Rxy_upper_right_corners, "Rxy_upper_right_corners");
  bout_mesh->get(Rxy_upper_left_corners, "Rxy_upper_left_corners");
  bout_mesh->get(Zxy_lower_left_corners, "Zxy_corners");
  bout_mesh->get(Zxy_lower_right_corners, "Zxy_lower_right_corners");
  bout_mesh->get(Zxy_upper_right_corners, "Zxy_upper_right_corners");
  bout_mesh->get(Zxy_upper_left_corners, "Zxy_upper_left_corners");
  Field2D ivertex_lower_left_corners;
  Field2D ivertex_lower_right_corners;
  Field2D ivertex_upper_right_corners;
  Field2D ivertex_upper_left_corners;
  if (!use_cxx_ivertex) {
    bout_mesh->get(ivertex_lower_left_corners, "ivertex_lower_left_corners");
    bout_mesh->get(ivertex_lower_right_corners, "ivertex_lower_right_corners");
    bout_mesh->get(ivertex_upper_right_corners, "ivertex_upper_right_corners");
    bout_mesh->get(ivertex_upper_left_corners, "ivertex_upper_left_corners");
  }
  // local number of x cells, excluding guards
  int Nx = bout_mesh->xend - bout_mesh->xstart + 1;
  // local number of y cells, excluding guards
  int Ny = bout_mesh->yend - bout_mesh->ystart + 1;
  // output << "Nx " + std::to_string(Nx) + "Ny " + std::to_string(Ny) << "\n";
  // output << "Got here -1 \n";

  PETSCCHK(PetscInitializeNoArguments());
  auto sycl_target = std::make_shared<SYCLTarget>(0, PETSC_COMM_WORLD);
  const int mpi_size = sycl_target->comm_pair.size_parent;
  const int mpi_rank = sycl_target->comm_pair.rank_parent;
  // output << "Got here 0 \n";
  // global number of physical nonunique vertices stored in hypnotoad datasets
  const int N_nonunique_vertices = mpi_size * Nx * Ny;
  // arrays to fill with local data
  std::vector<double> local_Z_lower_left_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> local_R_lower_left_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> local_Z_lower_right_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> local_R_lower_right_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> local_Z_upper_right_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> local_R_upper_right_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> local_Z_upper_left_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> local_R_upper_left_vertices(N_nonunique_vertices, 0.0);
  // arrays to receive the summed data across ranks
  std::vector<double> global_Z_lower_left_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> global_R_lower_left_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> global_Z_lower_right_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> global_R_lower_right_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> global_Z_upper_right_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> global_R_upper_right_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> global_Z_upper_left_vertices(N_nonunique_vertices, 0.0);
  std::vector<double> global_R_upper_left_vertices(N_nonunique_vertices, 0.0);
  // fill these vectors with vertex values from the local rank
  // at indices determined by the local rank
  int icxy = Nx * Ny * mpi_rank;
  for (int ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (int iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      local_R_lower_left_vertices.at(icxy) = Rxy_lower_left_corners(ix, iy);
      local_Z_lower_left_vertices.at(icxy) = Zxy_lower_left_corners(ix, iy);
      local_R_lower_right_vertices.at(icxy) = Rxy_lower_right_corners(ix, iy);
      local_Z_lower_right_vertices.at(icxy) = Zxy_lower_right_corners(ix, iy);
      local_R_upper_right_vertices.at(icxy) = Rxy_upper_right_corners(ix, iy);
      local_Z_upper_right_vertices.at(icxy) = Zxy_upper_right_corners(ix, iy);
      local_R_upper_left_vertices.at(icxy) = Rxy_upper_left_corners(ix, iy);
      local_Z_upper_left_vertices.at(icxy) = Zxy_upper_left_corners(ix, iy);
      icxy++;
    }
  }
  // Perform Allreduce (sum) to get knowledge of vertices to all ranks
  MPICHK(MPI_Allreduce(local_R_lower_left_vertices.data(),
                       global_R_lower_left_vertices.data(), N_nonunique_vertices,
                       MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
  MPICHK(MPI_Allreduce(local_Z_lower_left_vertices.data(),
                       global_Z_lower_left_vertices.data(), N_nonunique_vertices,
                       MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
  MPICHK(MPI_Allreduce(local_R_lower_right_vertices.data(),
                       global_R_lower_right_vertices.data(), N_nonunique_vertices,
                       MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
  MPICHK(MPI_Allreduce(local_Z_lower_right_vertices.data(),
                       global_Z_lower_right_vertices.data(), N_nonunique_vertices,
                       MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
  MPICHK(MPI_Allreduce(local_R_upper_right_vertices.data(),
                       global_R_upper_right_vertices.data(), N_nonunique_vertices,
                       MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
  MPICHK(MPI_Allreduce(local_Z_upper_right_vertices.data(),
                       global_Z_upper_right_vertices.data(), N_nonunique_vertices,
                       MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
  MPICHK(MPI_Allreduce(local_R_upper_left_vertices.data(),
                       global_R_upper_left_vertices.data(), N_nonunique_vertices,
                       MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
  MPICHK(MPI_Allreduce(local_Z_upper_left_vertices.data(),
                       global_Z_upper_left_vertices.data(), N_nonunique_vertices,
                       MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
  // if (mpi_rank == 0) {
  //     std::cout << "Result of Allreduce (sum): ";
  //     for (double val : global_R_lower_left_vertices) {
  //         std::cout << val << " ";
  //     }
  //     std::cout << std::endl;
  //     std::cout << "N_nonunique_vertices=" << N_nonunique_vertices << std::endl;
  // }
  // Now dynamically determine a list of unique vertex points
  // constant to give us a vector that can definitely contain all points in the global
  // lists
  const int N_global_nonunique_vertices = 4 * mpi_size * Nx * Ny;
  std::vector<double> global_Z_vertices_buffer(N_global_nonunique_vertices, 0.0);
  std::vector<double> global_R_vertices_buffer(N_global_nonunique_vertices, 0.0);
  // fill the buffer vectors, checking each time if the point is unique
  // first point, outside loop
  global_Z_vertices_buffer.at(0) = global_Z_lower_left_vertices.at(0);
  global_R_vertices_buffer.at(0) = global_R_lower_left_vertices.at(0);
  int N_unique = 1; // we have one unique point in the buffer
  const double zero = 1.0e-8;
  // loop over lower left vertices
  collect_unique_points(global_Z_vertices_buffer, global_R_vertices_buffer, N_unique,
                        zero, global_Z_lower_left_vertices, global_R_lower_left_vertices);
  // loop over lower right vertices
  collect_unique_points(global_Z_vertices_buffer, global_R_vertices_buffer, N_unique,
                        zero, global_Z_lower_right_vertices,
                        global_R_lower_right_vertices);
  // loop over upper right vertices
  collect_unique_points(global_Z_vertices_buffer, global_R_vertices_buffer, N_unique,
                        zero, global_Z_upper_right_vertices,
                        global_R_upper_right_vertices);
  // loop over upper left vertices
  collect_unique_points(global_Z_vertices_buffer, global_R_vertices_buffer, N_unique,
                        zero, global_Z_upper_left_vertices, global_R_upper_left_vertices);
  // now make a vector of the size N_unique and fill from the buffer
  std::vector<double> global_Z_vertices(N_unique, 0.0);
  std::vector<double> global_R_vertices(N_unique, 0.0);
  for (int iv = 0; iv < N_unique; iv++) {
    global_Z_vertices.at(iv) = global_Z_vertices_buffer.at(iv);
    global_R_vertices.at(iv) = global_R_vertices_buffer.at(iv);
  }
  if (mpi_rank == 0) {
    std::cout << "Result of vertex collection: ";
    // for (int iv=0; iv<N_unique; iv++) {
    //     std::cout << "(" << global_R_vertices.at(iv) << ", " <<
    //     global_Z_vertices.at(iv) << ") ";
    // }
    // std::cout << std::endl;
    std::cout << "N_unique=" << N_unique << std::endl;
  }
  // ivertex arrays made in cxx, initialise with -1 index
  Field2D ivertex_lower_left_corners_cxx{-1, bout_mesh};
  Field2D ivertex_lower_right_corners_cxx{-1, bout_mesh};
  Field2D ivertex_upper_right_corners_cxx{-1, bout_mesh};
  Field2D ivertex_upper_left_corners_cxx{-1, bout_mesh};
  // now fill ivertex_corners arrays
  RZ_to_ivertex_vector(ivertex_lower_left_corners_cxx, global_Z_vertices,
                       global_R_vertices, zero, bout_mesh, Rxy_lower_left_corners,
                       Zxy_lower_left_corners);
  RZ_to_ivertex_vector(ivertex_lower_right_corners_cxx, global_Z_vertices,
                       global_R_vertices, zero, bout_mesh, Rxy_lower_right_corners,
                       Zxy_lower_right_corners);
  RZ_to_ivertex_vector(ivertex_upper_right_corners_cxx, global_Z_vertices,
                       global_R_vertices, zero, bout_mesh, Rxy_upper_right_corners,
                       Zxy_upper_right_corners);
  RZ_to_ivertex_vector(ivertex_upper_left_corners_cxx, global_Z_vertices,
                       global_R_vertices, zero, bout_mesh, Rxy_upper_left_corners,
                       Zxy_upper_left_corners);

  // First we setup the topology of the mesh.
  PetscInt num_cells_owned = Nx * Ny;
  // std::vector<double> cells(4*num_cells_owned);
  std::vector<PetscInt> cells;
  if (use_cxx_ivertex) {
    cells_definition_from_RZ_ivertex(
        cells, bout_mesh, Rxy_lower_left_corners, Rxy_lower_right_corners,
        Rxy_upper_right_corners, Rxy_upper_left_corners, Zxy_lower_left_corners,
        Zxy_lower_right_corners, Zxy_upper_right_corners, Zxy_upper_left_corners,
        ivertex_lower_left_corners_cxx, ivertex_lower_right_corners_cxx,
        ivertex_upper_right_corners_cxx, ivertex_upper_left_corners_cxx);
  } else {
    cells = cells_definition_from_RZ_ivertex(
        cells, bout_mesh, Rxy_lower_left_corners, Rxy_lower_right_corners,
        Rxy_upper_right_corners, Rxy_upper_left_corners, Zxy_lower_left_corners,
        Zxy_lower_right_corners, Zxy_upper_right_corners, Zxy_upper_left_corners,
        ivertex_lower_left_corners, ivertex_lower_right_corners,
        ivertex_upper_right_corners, ivertex_upper_left_corners);
  }
  // create nvertices, global_vertex_list_R, global_vertex_list_z variables
  int nvertices;
  std::vector<double> global_vertex_list_R;
  std::vector<double> global_vertex_list_Z;
  if (use_cxx_ivertex) {
    nvertices = N_unique;
    global_vertex_list_R = global_R_vertices;
    global_vertex_list_Z = global_Z_vertices;
  } else {
    load_vertex_information_from_netcdf(nvertices, global_vertex_list_R,
                                        global_vertex_list_Z);
  }
  // number of vertices to keep per process for passing to
  // DMPlexCreateFromCellListParallelPetsc
  int nvertex_per_process = std::floor(nvertices / mpi_size);
  int nvertex_remainder = nvertices - mpi_size * nvertex_per_process;
  int nvertex_this_process = nvertex_per_process;
  // include the remaining vertices on the last rank
  if (mpi_rank == mpi_size - 1) {
    nvertex_this_process += nvertex_remainder;
  }
  // starting vertex index
  int ivertex_minimum = mpi_rank * nvertex_per_process;
  int ivertex_maximum = mpi_rank * nvertex_per_process + nvertex_this_process - 1;
  /*
   * Each rank owns a contiguous block of global indices. We label our indices
   * lexicographically (row-wise). Sorting out the global vertex indexing is
   * probably one of the more tedious parts.
   */
  PetscInt num_vertices_owned = nvertex_this_process;

  /*
   * Create the coordinates for the block of vertices we pass to petsc. For an
   * existing mesh in memory this step will probably involve some MPI
   * communication to gather the blocks of coordinates on the ranks which pass
   * them to PETSc.
   */
  std::vector<PetscScalar> vertex_coords(num_vertices_owned * 2);
  // shift due to differing rank
  int ishift;
  for (int iv = 0; iv < nvertex_this_process; iv++) {
    ishift = mpi_rank * nvertex_per_process;
    vertex_coords.at(iv * 2 + 0) = global_vertex_list_R[iv + ishift];
    vertex_coords.at(iv * 2 + 1) = global_vertex_list_Z[iv + ishift];
  }
  // This DM will contain the DMPlex after we call the creation routine.
  DM dm;
  // Create the DMPlex from the cells and coordinates.
  PETSCCHK(DMPlexCreateFromCellListParallelPetsc(
      PETSC_COMM_WORLD, 2, num_cells_owned, num_vertices_owned, PETSC_DECIDE, 4,
      PETSC_TRUE, cells.data(), 2, vertex_coords.data(), NULL, NULL, &dm));

  // Label all of the boundary faces with 100 in the "Face Sets" label by using
  // the helper function label_all_dmplex_boundaries.
  PetscInterface::label_all_dmplex_boundaries(dm, PetscInterface::face_sets_label, 100);

  // // Label subsections of the boundary by specifing pairs of vertices and using
  // // the label_dmplex_edges helper function.
  // std::vector<PetscInt> vertex_starts, vertex_ends, edge_labels;

  // if (mpi_rank == mpi_size - 1) {
  //   // Top edge
  //   for (int px = 0; px < mpi_size; px++) {
  //     const PetscInt tx = (mpi_size + 1) * mpi_size + px;
  //     vertex_starts.push_back(tx);
  //     vertex_ends.push_back(tx + 1);
  //     // Label the top edge with label 200
  //     edge_labels.push_back(200);
  //   }
  // }

  // PetscInterface::label_dmplex_edges(dm, PetscInterface::face_sets_label,
  //                                    vertex_starts, vertex_ends, edge_labels);

  // save a HDF5 file containing the DM for diagnostics
  PetscViewer viewer;
  // Set a name for the DMPlex object (important for HDF5)
  PetscObjectSetName((PetscObject)dm, dmplex_name.c_str());
  // Create an HDF5 viewer
  PetscViewerHDF5Open(PETSC_COMM_WORLD, dmplex_h5_filename.c_str(), FILE_MODE_WRITE,
                      &viewer);
  // Set viewer format to PETSC_VIEWER_HDF5_PETSC for compatibility
  PetscViewerPushFormat(viewer, PETSC_VIEWER_HDF5_PETSC);
  // Save the DMPlex to the HDF5 file
  DMView(dm, viewer);
  // Clean up
  PetscViewerDestroy(&viewer);
  output << "Finished DMPlex creation and diagnostic \n";
  output << "Begin particle push \n";
  // get data from BOUT for forces in particle push
  Field2D phi{bout_mesh};
  Field2D ex{bout_mesh};
  Field2D ey{bout_mesh};
  auto& opt = Options::root();
  phi = opt["mesh"]["phi"].as<Field2D>();
  ex = opt["mesh"]["ex"].as<Field2D>();
  ey = opt["mesh"]["ey"].as<Field2D>();

  /*
   *
   *
   *
   *
   *
   * Below here we setup basic advection on the dmplex
   *
   *
   *
   *
   *
   *
   */
  {
    const int ndim = 2;
    const int npart_per_cell =
        Options::root()["neso_particles"]["npart_per_cell"].withDefault(1);
    const REAL dt = Options::root()["neso_particles"]["dt"].withDefault(0.01);
    const int nsteps = Options::root()["neso_particles"]["nsteps"].withDefault(10);
    Field2D density{bout_mesh};
    bout_mesh->get(density, "density", 0.0, false);
    bout_mesh->communicate(density);

    // Create a mesh interface from the DM
    auto neso_mesh =
        std::make_shared<PetscInterface::DMPlexInterface>(dm, 0, MPI_COMM_WORLD);
    // Create a mapper for mapping particles into cells.
    auto mapper =
        std::make_shared<PetscInterface::DMPlexLocalMapper>(sycl_target, neso_mesh);
    // Create a domain from the neso_mesh and the mapper.
    auto domain = std::make_shared<Domain>(neso_mesh, mapper);

    // Create the particle properties (note that if you are using the Reactions
    // project it has its owne particle spec builder).
    ParticleSpec particle_spec{ParticleProp(Sym<REAL>("P"), ndim, true),
                               ParticleProp(Sym<REAL>("V"), ndim),
                               ParticleProp(Sym<REAL>("F"), ndim),
                               ParticleProp(Sym<REAL>("Q"), 1),
                               ParticleProp(Sym<REAL>("TSP"), 2),
                               ParticleProp(Sym<INT>("CELL_ID"), 1, true),
                               ParticleProp(Sym<INT>("ID"), 1)};

    // Create a Particle group with our specied particle properties.
    auto A = std::make_shared<ParticleGroup>(domain, particle_spec, sycl_target);

    // Create some particle data

    std::mt19937 rng_pos(52234234 + mpi_rank);
    std::mt19937 rng_vel(52234231 + mpi_rank);
    std::vector<std::vector<double>> positions;
    std::vector<int> cells;

    uniform_within_dmplex_cells(neso_mesh, npart_per_cell, positions, cells, &rng_pos);

    const int N_actual = cells.size();
    auto velocities =
        NESO::Particles::normal_distribution(N_actual, 2, 0.0, 1.0, rng_vel);

    int id_offset = 0;
    MPICHK(MPI_Exscan(&N_actual, &id_offset, 1, MPI_INT, MPI_SUM,
                      sycl_target->comm_pair.comm_parent));

    // This is host space to create particle data in before pushing the
    // particles into the ParticleGroup
    ParticleSet initial_distribution(N_actual, particle_spec);
    for (int px = 0; px < N_actual; px++) {
      for (int dimx = 0; dimx < ndim; dimx++) {
        initial_distribution[Sym<REAL>("P")][px][dimx] = positions[dimx][px];
        initial_distribution[Sym<REAL>("V")][px][dimx] = velocities[dimx][px];
      }
      // initial_distribution[Sym<REAL>("F")][px][0] = ex_cellorder(cells.at(px));
      // initial_distribution[Sym<REAL>("F")][px][1] = ey_cellorder(cells.at(px));
      initial_distribution[Sym<REAL>("F")][px][0] = 0.0;
      initial_distribution[Sym<REAL>("F")][px][1] = 0.0;
      initial_distribution[Sym<INT>("CELL_ID")][px][0] = cells.at(px);
      initial_distribution[Sym<INT>("ID")][px][0] = px + id_offset;
      initial_distribution[Sym<REAL>("Q")][px][0] = 1.0;
    }

    // Add the new particles to the particle group
    A->add_particles_local(initial_distribution);

    // make pointer to projection object
    auto dg0 = std::make_shared<PetscInterface::DMPlexProjectEvaluateDG>(
        neso_mesh, sycl_target, "DG", 0);
    // set F from a field from BOUT
    std::vector<REAL> h_project2(num_cells_owned * 2);
    // h_project2.reserve(num_cells_owned * 2);
    // this call should allocate h_project2
    // dg0->set_dofs(2, h_project2);
    PetscInt ixy = 0;
    for (PetscInt ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
      for (PetscInt iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
        h_project2.at(2 * ixy) = ex(ix, iy);
        h_project2.at(2 * ixy + 1) = ey(ix, iy);
        ixy++;
      }
    }
    // now copy the data to internal variables
    dg0->set_dofs(2, h_project2);
    // set the data from internal variables into F
    dg0->evaluate(A, Sym<REAL>("F"));

    // Create the boundary interaction objects
    std::map<PetscInt, std::vector<PetscInt>> boundary_groups;
    // boundary_groups[1] = {100, 200};
    boundary_groups[1] = {100};

    auto b2d = std::make_shared<PetscInterface::BoundaryInteraction2D>(
        sycl_target, neso_mesh, boundary_groups);
    auto reflection = std::make_shared<BoundaryReflection>(ndim, 1.0e-10);

    auto lambda_apply_boundary_conditions = [&](auto aa) {
      auto sub_groups = b2d->post_integration(aa);
      for (auto& gx : sub_groups) {
        reflection->execute(gx.second, Sym<REAL>("P"), Sym<REAL>("V"), Sym<REAL>("TSP"),
                            b2d->previous_position_sym);
      }
    };
    auto lambda_apply_timestep_reset = [&](auto aa) {
      particle_loop(
          aa,
          [=](auto TSP) {
            TSP.at(0) = 0.0;
            TSP.at(1) = 0.0;
          },
          Access::write(Sym<REAL>("TSP")))
          ->execute();
    };
    auto lambda_apply_advection_step =
        [=](ParticleSubGroupSharedPtr iteration_set) -> void {
      particle_loop(
          "euler_advection", iteration_set,
          [=](auto F, auto V, auto P, auto TSP) {
            const REAL dt_left = dt - TSP.at(0);
            if (dt_left > 0.0) {
              P.at(0) += dt_left * V.at(0);
              P.at(1) += dt_left * V.at(1);
              // update V from fixed F
              // n.b. a different advection scheme should be
              // used for more than 1st order accuracy
              V.at(0) += dt_left * F.at(0);
              V.at(1) += dt_left * F.at(1);
              TSP.at(0) = dt;
              TSP.at(1) = dt_left;
            }
          },
          Access::read(Sym<REAL>("F")), Access::write(Sym<REAL>("V")),
          Access::write(Sym<REAL>("P")), Access::write(Sym<REAL>("TSP")))
          ->execute();
    };
    auto lambda_pre_advection = [&](auto aa) { b2d->pre_integration(aa); };
    auto lambda_find_partial_moves = [&](auto aa) {
      return static_particle_sub_group(
          aa, [=](auto TSP) { return TSP.at(0) < dt; }, Access::read(Sym<REAL>("TSP")));
    };
    auto lambda_partial_moves_remaining = [&](auto aa) -> bool {
      const int size = get_npart_global(aa);
      ;
      return size > 0;
    };
    auto lambda_apply_timestep = [&](auto aa) {
      lambda_apply_timestep_reset(aa);
      lambda_pre_advection(aa);
      lambda_apply_advection_step(aa);
      lambda_apply_boundary_conditions(aa);
      aa = lambda_find_partial_moves(aa);
      while (lambda_partial_moves_remaining(aa)) {
        lambda_pre_advection(aa);
        lambda_apply_advection_step(aa);
        lambda_apply_boundary_conditions(aa);
        aa = lambda_find_partial_moves(aa);
      }
    };
    // for(int ix = bout_mesh->xstart; ix<= bout_mesh->xend; ix++){
    //   for(int iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++){
    //       //for(int iz=0; iz < bout_mesh->LocalNz; iz++){
    //         std::string string_count = std::string("(") + std::to_string(ix) +
    //         std::string(",") + std::to_string(iy) + std::string(")"); output <<
    //         string_count + std::string(": ") + std::to_string(phi(ix,iy)) +
    //         std::string("; ");
    //       //}
    //   }
    //   output << "\n";
    // }
    // uncomment to write a trajectory
    H5Part h5part("traj_reflection_dmplex_example.h5part", A, Sym<REAL>("P"),
                  Sym<REAL>("V"));

    // get a density by projecting the particle property Q to the bout_mesh
    dg0->project(A, Sym<REAL>("Q"));
    std::vector<REAL> h_project1;
    dg0->get_dofs(1, h_project1);
    PetscInt ic = 0;
    for (PetscInt ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
      for (PetscInt iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
        density(ix, iy) = h_project1.at(ic);
        ic++;
      }
    }
    // this fills internal guards
    bout_mesh->communicate(density);
    // apply boundary conditions to fill external guards
    // density.applyBoundary();
    // extrapolate -> Neumann
    // BOUT++ may need to do something with density Field2D guards cells
    // print density to screen to show non-trivial result
    // compare to 1/J*dx*dy*dz -> at the initial time we have 1 particle per cell
    // so the density is 1/Cell_volume
    // for(int ix = bout_mesh->xstart; ix<= bout_mesh->xend; ix++){
    //   for(int iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++){
    //       //for(int iz=0; iz < bout_mesh->LocalNz; iz++){
    //         std::string string_count = std::string("(") + std::to_string(ix) +
    //         std::string(",") + std::to_string(iy) + std::string(")"); output <<
    //         string_count + std::string(": ") + std::to_string(density(ix,iy)) +
    //         std::string("; ") +
    //         std::to_string(1.0/(coord->J(ix,iy)*coord->dx(ix,iy)*coord->dy(ix,iy)*coord->dz(ix,iy)));
    //       //}
    //   }
    //   output << "\n";
    // }

    for (int stepx = 0; stepx < nsteps; stepx++) {
      // nprint("step:", stepx);
      output << "step:" << std::to_string(stepx) << std::endl;
      A->hybrid_move();
      A->cell_move();
      lambda_apply_timestep(static_particle_sub_group(A));

      // uncomment to write a trajectory
      h5part.write();
    }

    // uncomment to write a trajectory
    h5part.close();

    // Boundary interaction objects require a free call.
    b2d->free();
    // NESO-Particles neso_mesh objects must have free called on them.
    neso_mesh->free();
  }
  PETSCCHK(DMDestroy(&dm));
  PETSCCHK(PetscFinalize());
  sycl_target->free();

  if (MPI_Finalize() != MPI_SUCCESS) {
    std::cout << "ERROR: MPI_Finalize != MPI_SUCCESS" << std::endl;
    return -1;
  }
  return 0;
}

#endif
