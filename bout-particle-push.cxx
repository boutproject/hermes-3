#include "bout/bout.hxx"
#include "bout/field2d.hxx"
#include "bout/output.hxx"
#include "bout/petsclib.hxx"
#include <bout/field_factory.hxx>
#include <fmt/core.h>
#include <iostream>
#include <memory>
#include <neso_particles.hpp>
#include <neso_particles/compute_target.hpp>
#include <neso_particles/external_interfaces/petsc/petsc_interface.hpp>
#include <neso_particles/typedefs.hpp>
#include <netcdf>
#include <petscsystypes.h>
#include <petscviewerhdf5.h>
#include <string>
#include <vector>
// for reactions integration
#include <reactions/reactions.hpp>

#ifndef NESO_PARTICLES_PETSC
static_assert(false, "NESO-Particles was installed without PETSc support.");
#else

using namespace NESO::Particles;
using namespace VANTAGE::Reactions;

template <typename T, typename U>
inline void ASSERT_EQ(T t, U u) {
  NESOASSERT(t == u, "A check failed.");
}

void collect_unique_points(std::vector<double>& global_Z_vertices_buffer,
                           std::vector<double>& global_R_vertices_buffer, int& N_unique,
                           const double& tolerance,
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
              < tolerance
          && abs(global_R_hypnotoad_vertices.at(iv)
                 - global_R_vertices_buffer.at(iunique))
                 < tolerance) {
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
                          std::vector<double>& global_R_vertices, const double& tolerance,
                          Mesh*& bout_mesh, Field2D& Rxy_corners, Field2D& Zxy_corners) {
  int Nvertex = global_Z_vertices.size();
  bool index_found;
  // loop over points that are not guard cells
  for (int ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (int iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      index_found = false;
      for (int iv = 0; iv < Nvertex; iv++) {
        if (abs(global_Z_vertices.at(iv) - Zxy_corners(ix, iy)) < tolerance
            && abs(global_R_vertices.at(iv) - Rxy_corners(ix, iy)) < tolerance) {
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

DM create_dmplex_from_Bout_mesh(Mesh* bout_mesh, std::shared_ptr<SYCLTarget> sycl_target){
    bool use_cxx_ivertex = Options::root()["mesh"]["use_cxx_ivertex"].withDefault(false);
    std::string dmplex_name =
        Options::root()["mesh"]["dmplex_name"].withDefault("hypnotoad_dmplex_mesh");
    std::string dmplex_h5_filename =
        Options::root()["mesh"]["dmplex_h5_filename"].withDefault(
            "hypnotoad_dmplex_mesh_output.h5");
    output << fmt::format("Using option use_cxx_ivertex = {}", use_cxx_ivertex)
          << std::endl;
    bout_mesh->load();
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

    // PETSCCHK(PetscInitializeNoArguments());
    // auto sycl_target = std::make_shared<SYCLTarget>(0, PETSC_COMM_WORLD);
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
    const double tolerance = 1.0e-8;
    // loop over lower left vertices
    collect_unique_points(global_Z_vertices_buffer, global_R_vertices_buffer, N_unique,
                          tolerance, global_Z_lower_left_vertices, global_R_lower_left_vertices);
    // loop over lower right vertices
    collect_unique_points(global_Z_vertices_buffer, global_R_vertices_buffer, N_unique,
                          tolerance, global_Z_lower_right_vertices,
                          global_R_lower_right_vertices);
    // loop over upper right vertices
    collect_unique_points(global_Z_vertices_buffer, global_R_vertices_buffer, N_unique,
                          tolerance, global_Z_upper_right_vertices,
                          global_R_upper_right_vertices);
    // loop over upper left vertices
    collect_unique_points(global_Z_vertices_buffer, global_R_vertices_buffer, N_unique,
                          tolerance, global_Z_upper_left_vertices, global_R_upper_left_vertices);
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
                        global_R_vertices, tolerance, bout_mesh, Rxy_lower_left_corners,
                        Zxy_lower_left_corners);
    RZ_to_ivertex_vector(ivertex_lower_right_corners_cxx, global_Z_vertices,
                        global_R_vertices, tolerance, bout_mesh, Rxy_lower_right_corners,
                        Zxy_lower_right_corners);
    RZ_to_ivertex_vector(ivertex_upper_right_corners_cxx, global_Z_vertices,
                        global_R_vertices, tolerance, bout_mesh, Rxy_upper_right_corners,
                        Zxy_upper_right_corners);
    RZ_to_ivertex_vector(ivertex_upper_left_corners_cxx, global_Z_vertices,
                        global_R_vertices, tolerance, bout_mesh, Rxy_upper_left_corners,
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
    return dm;
}

void calculate_density_in_place(Field2D& density,
      std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG>& dg0,
      std::shared_ptr<ParticleGroup>& A_particle_group,
      std::vector<double>& h_project1) {
  Mesh* bout_mesh = density.getMesh();
  // get a density by projecting the particle property WEIGHT to the bout_mesh
  dg0->project(A_particle_group, Sym<REAL>("WEIGHT"));
  // std::vector<REAL> h_project1;
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
}

double calculate_total_mass(Field2D& density, std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh){
  double local_mass=0.0;
  double total_mass=0.0;
  Mesh* bout_mesh = density.getMesh();
  PetscInt ic = 0;
  for (PetscInt ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (PetscInt iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      local_mass += density(ix, iy)*neso_mesh->dmh->get_cell_volume(ic);
      ic++;
    }
  }
  MPICHK(MPI_Allreduce(&local_mass,
                       &total_mass, 1,
                        MPI_DOUBLE, MPI_SUM, BoutComm::get()));
  return total_mass;
}

Options initialise_diagnostics(Field2D& density,
      std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG>& dg0,
      std::shared_ptr<ParticleGroup>& A_particle_group,
      std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
      std::vector<double>& h_project1,
      std::string particle_data_filename){
  // Options object to use to write out diagnostic data of fluid quantities
  Options bout_output_data;
  calculate_density_in_place(density, dg0, A_particle_group, h_project1);
  bout_output_data["neutral_density"] = density;
  // Set the time attribute
  bout_output_data["neutral_density"].attributes["time_dimension"] = "t";
  bout_output_data["total_mass"] = calculate_total_mass(density, neso_mesh);
  bout_output_data["total_mass"].attributes["time_dimension"] = "t";
  // std::string particle_data_filename = fmt::format("bout_particle_moments_{}.nc",mpi_rank);
  bout::OptionsIO::create(particle_data_filename)->write(bout_output_data);
  return bout_output_data;
}

void update_diagnostics(Field2D& density,
      std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG>& dg0,
      std::shared_ptr<ParticleGroup>& A_particle_group,
      std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
      std::vector<double>& h_project1,
      Options& bout_output_data,
      std::string particle_data_filename){
  // update density and write
  calculate_density_in_place(density, dg0, A_particle_group, h_project1);
  bout_output_data["neutral_density"] = density;
  bout_output_data["total_mass"] = calculate_total_mass(density, neso_mesh);
  // Append data to file
  bout::OptionsIO::create({{"file", particle_data_filename}, {"append", true}})->write(bout_output_data);
}

void set_initial_particle_weights(Field2D& initial_density,
      std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG>& dg0,
      std::shared_ptr<ParticleGroup>& A_particle_group,
      std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
      std::vector<double>& h_project1){
  Mesh* bout_mesh = initial_density.getMesh();
  PetscInt ixy = 0;
  for (PetscInt ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (PetscInt iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      // particle_weights are copied to all particles in this cell so
      // we multiply the initial density by the volume to get particle number,
      // then divide by the number of marker particles per cell
      const REAL cell_volume = neso_mesh->dmh->get_cell_volume(ixy);
      const INT nmarkers_per_cell = A_particle_group->get_npart_cell(ixy);
      const REAL particle_weights = initial_density(ix, iy)*cell_volume/nmarkers_per_cell;
      h_project1.at(ixy) = particle_weights;
      ixy++;
    }
  }
  // now copy the data to internal variables
  dg0->set_dofs(1, h_project1);
  // set the data from internal variables into the weights
  dg0->evaluate(A_particle_group, Sym<REAL>("WEIGHT"));
}

void check_cell_volumes(
      std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
      Mesh*& bout_mesh){
  Coordinates* coord = bout_mesh->getCoordinates();
  PetscInt ixy = 0;
  const REAL tolerance = 1.0e-12;
  for (PetscInt ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (PetscInt iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      const REAL bout_cell_volume = coord->J(ix,iy)*
                                coord->dx(ix,iy)*coord->dy(ix,iy);
      const REAL neso_cell_volume = neso_mesh->dmh->get_cell_volume(ixy);
      const bool volumes_match = (abs(bout_cell_volume - neso_cell_volume) < tolerance);
      // exit if we fail to find a match
      NESOASSERT(volumes_match,
                 fmt::format("BOUT++ mesh volume {} does not match NESO-Particles mesh volume {} for ix = {} iy = {}", bout_cell_volume, neso_cell_volume, ix, iy));
      ixy++;
    }
  }
}

int main(int argc, char** argv) {
  // initialise_mpi(&argc, &argv);
  // attempt to call BOUT to
  // get information from a BOUT
  // mesh object
  // N.B. Comment the next three lines
  // to permit compilation as is
  BoutInitialise(argc, argv);
  Mesh* bout_mesh = Mesh::create(&Options::root()["mesh"]);
  PETSCCHK(PetscInitializeNoArguments());
  auto sycl_target = std::make_shared<SYCLTarget>(0, PETSC_COMM_WORLD);
  DM dm = create_dmplex_from_Bout_mesh(bout_mesh, sycl_target);
  output << "Begin particle push \n";
  // get data from BOUT.inp to assign particle weights as a fn of x,y
  auto& opt = Options::root();
  Field2D initial_density{bout_mesh};
  initial_density = opt["mesh"]["initial_density"].as<Field2D>();
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
    Field2D density = Field2D(0.0, bout_mesh);
    // Create a mesh interface from the DM
    auto neso_mesh =
        std::make_shared<PetscInterface::DMPlexInterface>(dm, 0, MPI_COMM_WORLD);
    // Create a mapper for mapping particles into cells.
    auto mapper =
        std::make_shared<PetscInterface::DMPlexLocalMapper>(sycl_target, neso_mesh);
    // Create a domain from the neso_mesh and the mapper.
    auto domain = std::make_shared<Domain>(neso_mesh, mapper);
    // Get the number of cells in the mesh owned on this process
    int num_cells_owned = neso_mesh->get_cell_count();
    // if requested, check that neso_mesh cell volumes are identical
    // to bout_mesh cell volumes, otherwise, exit.
    if (Options::root()["neso_particles"]["test_cell_volumes"].withDefault(false)){
      check_cell_volumes(neso_mesh, bout_mesh);
    }
    // create a Reactions particle spec
    auto particle_spec_builder = ParticleSpecBuilder(ndim);
    auto electron_species = Species("ELECTRON");
    auto main_species = Species("ION", 1.0, 0.0, 0);
    std::vector<Species> fluid_species = {electron_species, main_species};
    particle_spec_builder.add_particle_prop(Properties<REAL>(
        fluid_species, std::vector<int>{default_properties.temperature,
                                        default_properties.density,
                                        default_properties.source_energy,
                                        default_properties.source_density}));
    particle_spec_builder.add_particle_prop(
        Properties<REAL>(fluid_species,
                        std::vector<int>{default_properties.source_momentum}),
        ndim);
    ParticleSpec additional_props{
      ParticleProp(Sym<REAL>("TSP"), 2)};
    particle_spec_builder.add_particle_spec(additional_props);
    ParticleSpec particle_spec = particle_spec_builder.get_particle_spec();

    // Create a Particle group with our specied particle properties.
    auto A_particle_group = std::make_shared<ParticleGroup>(domain, particle_spec, sycl_target);

    // Create some particle data
    const int mpi_rank = sycl_target->comm_pair.rank_parent;
    std::mt19937 rng_pos(52234234 + mpi_rank);
    std::mt19937 rng_vel(52234231 + mpi_rank);
    std::vector<std::vector<double>> positions;
    std::vector<int> particle_cell_ids;

    uniform_within_dmplex_cells(neso_mesh, npart_per_cell, positions, particle_cell_ids, &rng_pos);

    const int N_actual = particle_cell_ids.size();
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
        initial_distribution[Sym<REAL>("POSITION")][px][dimx] = positions[dimx][px];
        initial_distribution[Sym<REAL>("VELOCITY")][px][dimx] = velocities[dimx][px];
      }
      initial_distribution[Sym<REAL>("ION_DENSITY")][px][0] = 1.0;
      initial_distribution[Sym<REAL>("ION_SOURCE_DENSITY")][px][0] = 1.0;
      initial_distribution[Sym<REAL>("ION_SOURCE_ENERGY")][px][0] = 0.5;
      initial_distribution[Sym<REAL>("ELECTRON_DENSITY")][px][0] = 1.0;
      initial_distribution[Sym<REAL>("ELECTRON_SOURCE_DENSITY")][px][0] = 1.0;
      initial_distribution[Sym<REAL>("ELECTRON_SOURCE_ENERGY")][px][0] = 0.5;
      for (int dimx = 0; dimx < ndim; dimx++) {
        initial_distribution[Sym<REAL>("ION_SOURCE_MOMENTUM")][px][dimx] = 0.0;
        initial_distribution[Sym<REAL>("ELECTRON_SOURCE_MOMENTUM")][px][dimx] = 0.0;
      }
      initial_distribution[Sym<INT>("CELL_ID")][px][0] = particle_cell_ids.at(px);
      initial_distribution[Sym<INT>("ID")][px][0] = px + id_offset;
      initial_distribution[Sym<REAL>("WEIGHT")][px][0] = 1.0;
    }
    // Add the new particles to the particle group
    A_particle_group->add_particles_local(initial_distribution);
    // make pointer to projection object
    auto dg0 = std::make_shared<PetscInterface::DMPlexProjectEvaluateDG>(
        neso_mesh, sycl_target, "DG", 0);
    const REAL iz_rate = Options::root()["VANTAGE_reactions"]["iz_rate"].withDefault(1.0);
    auto iz_rate_data = FixedRateData(iz_rate);
    main_species.set_id(0);
    auto ionisation_reaction = ElectronImpactIonisation<FixedRateData, FixedRateData>(
        A_particle_group->sycl_target, iz_rate_data, iz_rate_data, main_species,
        electron_species);

    // Reaction controllers
    const REAL remove_threshold = Options::root()["VANTAGE_reactions"]["remove_threshold"].withDefault(1.0e-10);
    const REAL merge_threshold = Options::root()["VANTAGE_reactions"]["merge_threshold"].withDefault(1.0e-2);

    auto remove_wrapper = std::make_shared<TransformationWrapper>(
      std::vector<std::shared_ptr<MarkingStrategy>>{
          make_marking_strategy<ComparisonMarkerSingle<REAL, LessThanComp>>(
              Sym<REAL>("WEIGHT"), remove_threshold)},
      std::dynamic_pointer_cast<TransformationStrategy>(std::make_shared<SimpleRemovalTransformationStrategy>()));

    auto merge_wrapper = std::make_shared<TransformationWrapper>(
      std::vector<std::shared_ptr<MarkingStrategy>>{
          make_marking_strategy<ComparisonMarkerSingle<REAL, LessThanComp>>(
              Sym<REAL>("WEIGHT"), merge_threshold)},
      make_transformation_strategy<MergeTransformationStrategy<ndim>>());

    auto child_transforms = std::vector{merge_wrapper,remove_wrapper};
    std::vector<std::shared_ptr<TransformationWrapper>> parent_transforms = child_transforms;

    auto reaction_controller =
        ReactionController(parent_transforms, child_transforms);
    // add ionisation to the controller
    reaction_controller.add_reaction(std::make_shared<decltype(ionisation_reaction)>(ionisation_reaction));

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
        reflection->execute(gx.second, Sym<REAL>("POSITION"), Sym<REAL>("VELOCITY"), Sym<REAL>("TSP"),
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
          [=](auto VELOCITY, auto POSITION, auto TSP) {
            const REAL dt_left = dt - TSP.at(0);
            if (dt_left > 0.0) {
              POSITION.at(0) += dt_left * VELOCITY.at(0);
              POSITION.at(1) += dt_left * VELOCITY.at(1);
              TSP.at(0) = dt;
              TSP.at(1) = dt_left;
            }
          },
          Access::read(Sym<REAL>("VELOCITY")),
          Access::write(Sym<REAL>("POSITION")),
          Access::write(Sym<REAL>("TSP")))
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
    // uncomment to write a trajectory
    H5Part h5part("traj_reflection_dmplex_example.h5part", A_particle_group, Sym<REAL>("POSITION"),
                  Sym<REAL>("VELOCITY"));

    // allocate buffer vector for scalar projection/evaluation of NESO-Particles properties
    std::vector<REAL> h_project1(num_cells_owned);
    // set weights from a Field2D from BOUT
    set_initial_particle_weights(initial_density, dg0, A_particle_group, neso_mesh, h_project1);
    // diagnose the initial condition
    std::string particle_data_filename = fmt::format("bout_particle_moments_{}.nc",mpi_rank);
    Options bout_output_data = initialise_diagnostics(density, dg0,A_particle_group,
        neso_mesh,h_project1, particle_data_filename);
    // begin timestepping
    for (int stepx = 0; stepx < nsteps; stepx++) {
      // nprint("step:", stepx);
      output << "step:" << std::to_string(stepx) << std::endl;
      A_particle_group->hybrid_move();
      A_particle_group->cell_move();
      lambda_apply_timestep(static_particle_sub_group(A_particle_group));
      // apply reactions
      reaction_controller.apply(A_particle_group, dt, ControllerMode::standard_mode);
      // uncomment to write a trajectory
      h5part.write();
      // uncomment to print particle info
      // A_particle_group->print(Sym<REAL>("POSITION"), Sym<INT>("ID"), Sym<REAL>("WEIGHT"));
      // diagnose timestep stepx
      update_diagnostics(density, dg0, A_particle_group, neso_mesh,
        h_project1, bout_output_data, particle_data_filename);
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
