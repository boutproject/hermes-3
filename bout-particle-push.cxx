#include "bout/field2d.hxx"
#include "bout/output.hxx"
#include "bout/petsclib.hxx"
#include <neso_particles.hpp>
#include <neso_particles/external_interfaces/petsc/petsc_interface.hpp>
#include <petscsystypes.h>
#include <string>
#include "bout/bout.hxx"
#include <bout/field_factory.hxx>
#include <netcdf>
#include <iostream>
#include <petscviewerhdf5.h>

#ifndef NESO_PARTICLES_PETSC
static_assert(false, "NESO-Particles was installed without PETSc support.");
#else

template<typename T, typename U>
inline void ASSERT_EQ(T t, U u){
  NESOASSERT(t==u, "A check failed.");
}


using namespace NESO::Particles;

int main(int argc, char **argv) {
  // initialise_mpi(&argc, &argv);
  // attempt to call BOUT to
  // get information from a BOUT
  // mesh object
  // N.B. Comment the next three lines
  // to permit compilation as is
  BoutInitialise(argc, argv);
  Mesh* mesh = Mesh::create(&Options::root()["mesh"]);
  mesh->load();
  Field2D Rxy_lower_left_corners;
  Field2D Rxy_lower_right_corners;
  Field2D Rxy_upper_right_corners;
  Field2D Rxy_upper_left_corners;
  Field2D Zxy_lower_left_corners;
  Field2D Zxy_lower_right_corners;
  Field2D Zxy_upper_right_corners;
  Field2D Zxy_upper_left_corners;
  Field2D ivertex_lower_left_corners;
  Field2D ivertex_lower_right_corners;
  Field2D ivertex_upper_right_corners;
  Field2D ivertex_upper_left_corners;
  //mesh->get(ivertex, "ivertex_lower_left_corners");
  mesh->get(Rxy_lower_left_corners, "Rxy_corners");
  mesh->get(Rxy_lower_right_corners, "Rxy_lower_right_corners");
  mesh->get(Rxy_upper_right_corners, "Rxy_upper_right_corners");
  mesh->get(Rxy_upper_left_corners, "Rxy_upper_left_corners");
  mesh->get(Zxy_lower_left_corners, "Zxy_corners");
  mesh->get(Zxy_lower_right_corners, "Zxy_lower_right_corners");
  mesh->get(Zxy_upper_right_corners, "Zxy_upper_right_corners");
  mesh->get(Zxy_upper_left_corners, "Zxy_upper_left_corners");
  mesh->get(ivertex_lower_left_corners, "ivertex_lower_left_corners");
  mesh->get(ivertex_lower_right_corners, "ivertex_lower_right_corners");
  mesh->get(ivertex_upper_right_corners, "ivertex_upper_right_corners");
  mesh->get(ivertex_upper_left_corners, "ivertex_upper_left_corners");
  // for(int ix = mesh->xstart; ix<= mesh->xend; ix++){
  //  for(int iy = mesh->ystart; iy <= mesh->yend; iy++){
  //      //for(int iz=0; iz < mesh->LocalNz; iz++){

  //        std::string string_count = std::string("(") + std::to_string(ix) + std::string(",") + std::to_string(iy) + std::string(")");
  //        output << string_count + std::string(": ") + std::to_string(static_cast<int>(ivertex_lower_left_corners(ix,iy))) + std::string("; ");
  //      //}
  //  }
  // output << "\n";
  // }
  // local number of x cells, excluding guards
  int Nx = mesh->xend - mesh->xstart + 1;
  // local number of y cells, excluding guards
  int Ny = mesh->yend - mesh->ystart + 1;
  // output << "Nx " + std::to_string(Nx) + "Ny " + std::to_string(Ny) << "\n";
  // output << "Got here -1 \n";

  PETSCCHK(PetscInitializeNoArguments());
  auto sycl_target = std::make_shared<SYCLTarget>(0, PETSC_COMM_WORLD);
  const int mpi_size = sycl_target->comm_pair.size_parent;
  const int mpi_rank = sycl_target->comm_pair.rank_parent;
  // output << "Got here 0 \n";

  // First we setup the topology of the mesh.
  std::vector<PetscReal> Z_vertices(4);
  std::vector<PetscReal> R_vertices(4);
  std::vector<PetscReal> theta_vertices(4);
  std::vector<PetscInt> i_vertices(4);
  std::vector<PetscInt> sort_indices(4);
  PetscReal ZZ;
  PetscReal ZZmid;
  PetscReal RR;
  PetscReal RRmid;
  PetscInt num_cells_owned = Nx*Ny;
  std::vector<PetscInt> cells;
  cells.reserve(num_cells_owned * 4);
  // output << "Got here 1 \n";
  // We are careful to list the vertices in counter clock-wise order, this might
  // matter.
  for (PetscInt ix = mesh->xstart; ix<= mesh->xend; ix++) {
    for (PetscInt iy = mesh->ystart; iy <= mesh->yend; iy++) {
      // collect data from Hypnotoad arrays
      i_vertices[0] =  static_cast<int>(ivertex_lower_left_corners(ix,iy));
      i_vertices[1] =  static_cast<int>(ivertex_lower_right_corners(ix,iy));
      i_vertices[2] =  static_cast<int>(ivertex_upper_right_corners(ix,iy));
      i_vertices[3] =  static_cast<int>(ivertex_upper_left_corners(ix,iy));

      R_vertices[0] =  Rxy_lower_left_corners(ix,iy);
      R_vertices[1] =  Rxy_lower_right_corners(ix,iy);
      R_vertices[2] =  Rxy_upper_right_corners(ix,iy);
      R_vertices[3] =  Rxy_upper_left_corners(ix,iy);

      Z_vertices[0] =  Zxy_lower_left_corners(ix,iy);
      Z_vertices[1] =  Zxy_lower_right_corners(ix,iy);
      Z_vertices[2] =  Zxy_upper_right_corners(ix,iy);
      Z_vertices[3] =  Zxy_upper_left_corners(ix,iy);
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
        theta_vertices[iv] = std::atan2(ZZ,RR);
      }
      // get the indices that sort the vertices in ascending order of theta
      for (PetscInt iv = 0; iv < 4; ++iv) {
        sort_indices[iv] = iv;
      }
      // output << "Before sort \n";
      // for (PetscInt iv = 0; iv < 4; ++iv) {
      //   output << std::to_string(i_vertices[iv]) << " " << std::to_string(theta_vertices[iv]) << " " << std::to_string(R_vertices[iv]) << " " << std::to_string(Z_vertices[iv]) << "\n";
      // }
      std::sort(sort_indices.begin(), sort_indices.end(),
                [&theta_vertices](int i, int j)
                  {
                      return theta_vertices[i] < theta_vertices[j];
                  });
      // output << "After sort \n";
      // for (PetscInt iv = 0; iv < 4; ++iv) {
      //   output << std::to_string(i_vertices[sort_indices[iv]]) << " " << std::to_string(theta_vertices[sort_indices[iv]]) << "\n";
      // }
      // fill cells using the sorted indices
      for (PetscInt iv = 0; iv < 4; ++iv) {
        cells.push_back(i_vertices[sort_indices[iv]]);
      }
      // // These are global indices not local indices.
      // PetscInt vertex_sw = cx + mpi_rank * (mpi_size + 1);
      // PetscInt vertex_se = vertex_sw + 1;
      // PetscInt vertex_ne = vertex_se + mpi_size + 1;
      // PetscInt vertex_nw = vertex_ne - 1;
      // cells.push_back(vertex_sw);
      // cells.push_back(vertex_se);
      // cells.push_back(vertex_ne);
      // cells.push_back(vertex_nw);
    }
  }

  // read data from netcdf for global vertices in mesh
  // Open the NetCDF file in read-only mode
  const std::string filename = "particle-push-data/expected_nonorthogonal.grd.vertex.nc";
  netCDF::NcFile dataFile(filename, netCDF::NcFile::read);

  // Get the variable
  std::string varName = "global_vertex_list_R";
  netCDF::NcVar dataVar = dataFile.getVar(varName);
  if (dataVar.isNull()) {
      std::cerr << "Variable '" << varName << "' not found in file." << std::endl;
      return -1;
  }
  std::vector<netCDF::NcDim> dims = dataVar.getDims();
  size_t nvertices = dims[0].getSize();

  // Read the data into a vector
  std::vector<double> global_vertex_list_R(nvertices);
  dataVar.getVar(global_vertex_list_R.data());

  // Get the variable
  varName = "global_vertex_list_Z";
  dataVar = dataFile.getVar(varName);
  if (dataVar.isNull()) {
      std::cerr << "Variable '" << varName << "' not found in file." << std::endl;
      return -1;
  }
  // Read the data into a vector
  std::vector<double> global_vertex_list_Z(nvertices);
  dataVar.getVar(global_vertex_list_Z.data());

  // // Print the data
  // std::cout << "Data from variable '" << varName << "':" << std::endl;
  // for (size_t i = 0; i < length; ++i) {
  //     std::cout << data[i] << " ";
  // }
  // std::cout << std::endl;
  // std::cout << "nvertices " << std::to_string(nvertices) << std::endl;
  // number of vertices to keep per process
  int nvertex_per_process = std::floor(nvertices/mpi_size);
  int nvertex_remainder = nvertices - mpi_size*nvertex_per_process;
  int nvertex_this_process = nvertex_per_process;
  // include the remaining vertices on the last rank
  if (mpi_rank == mpi_size - 1) {
    nvertex_this_process += nvertex_remainder;
  }
  // starting vertex index
  int ivertex_minimum = mpi_rank*nvertex_per_process;
  int ivertex_maximum = mpi_rank*nvertex_per_process + nvertex_this_process - 1;
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
  // output << "Got here 2 \n";
  // shift due to differing rank
  int ishift;
  for (int iv = 0; iv < nvertex_this_process; iv++) {
    ishift = mpi_rank*nvertex_per_process;
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
  PetscInterface::label_all_dmplex_boundaries(
      dm, PetscInterface::face_sets_label, 100);

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
  PetscObjectSetName((PetscObject)dm, "hypnotoad_dmplex_mesh");
  // Create an HDF5 viewer
  PetscViewerHDF5Open(PETSC_COMM_WORLD, "hypnotoad_dmplex_mesh_output.h5", FILE_MODE_WRITE, &viewer);
  // Set viewer format to PETSC_VIEWER_HDF5_PETSC for compatibility
  PetscViewerPushFormat(viewer, PETSC_VIEWER_HDF5_PETSC);
  // Save the DMPlex to the HDF5 file
  DMView(dm, viewer);
  // Clean up
  PetscViewerDestroy(&viewer);
  output << "Finished DMPlex creation and diagnostic \n";
  output << "Begin particle push \n";
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
    const int npart_per_cell = 1;
    const REAL dt = 0.1;
    const int nsteps = 3;

    // Create a mesh interface from the DM
    auto neso_mesh = std::make_shared<PetscInterface::DMPlexInterface>(
        dm, 0, MPI_COMM_WORLD);
    // Create a mapper for mapping particles into cells.
    auto mapper =
        std::make_shared<PetscInterface::DMPlexLocalMapper>(sycl_target, neso_mesh);
    // Create a domain from the neso_mesh and the mapper.
    auto domain = std::make_shared<Domain>(neso_mesh, mapper);

    // Create the particle properties (note that if you are using the Reactions
    // project it has its owne particle spec builder).
    ParticleSpec particle_spec{ParticleProp(Sym<REAL>("P"), ndim, true),
                               ParticleProp(Sym<REAL>("V"), ndim),
                               ParticleProp(Sym<REAL>("TSP"), 2),
                               ParticleProp(Sym<INT>("CELL_ID"), 1, true),
                               ParticleProp(Sym<INT>("ID"), 1)};

    // Create a Particle group with our specied particle properties.
    auto A =
        std::make_shared<ParticleGroup>(domain, particle_spec, sycl_target);

    // Create some particle data

    std::mt19937 rng_pos(52234234 + mpi_rank);
    std::mt19937 rng_vel(52234231 + mpi_rank);
    std::vector<std::vector<double>> positions;
    std::vector<int> cells;

    uniform_within_dmplex_cells(neso_mesh, npart_per_cell, positions, cells,
                                &rng_pos);

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
      initial_distribution[Sym<INT>("CELL_ID")][px][0] = cells.at(px);
      initial_distribution[Sym<INT>("ID")][px][0] = px + id_offset;
    }

    // Add the new particles to the particle group
    A->add_particles_local(initial_distribution);

    // Create the boundary interaction objects
    std::map<PetscInt, std::vector<PetscInt>> boundary_groups;
    boundary_groups[1] = {100, 200};

    auto b2d = std::make_shared<PetscInterface::BoundaryInteraction2D>(
        sycl_target, neso_mesh, boundary_groups);
    auto reflection = std::make_shared<BoundaryReflection>(ndim, 1.0e-10);

    auto lambda_apply_boundary_conditions = [&](auto aa) {
      auto sub_groups = b2d->post_integration(aa);
      for (auto &gx : sub_groups) {
        reflection->execute(gx.second, Sym<REAL>("P"), Sym<REAL>("V"),
                            Sym<REAL>("TSP"), b2d->previous_position_sym);
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
          [=](auto V, auto P, auto TSP) {
            const REAL dt_left = dt - TSP.at(0);
            if (dt_left > 0.0) {
              P.at(0) += dt_left * V.at(0);
              P.at(1) += dt_left * V.at(1);
              TSP.at(0) = dt;
              TSP.at(1) = dt_left;
            }
          },
          Access::read(Sym<REAL>("V")), Access::write(Sym<REAL>("P")),
          Access::write(Sym<REAL>("TSP")))
          ->execute();
    };
    auto lambda_pre_advection = [&](auto aa) { b2d->pre_integration(aa); };
    auto lambda_find_partial_moves = [&](auto aa) {
      return static_particle_sub_group(
          aa, [=](auto TSP) { return TSP.at(0) < dt; },
          Access::read(Sym<REAL>("TSP")));
    };
    auto lambda_partial_moves_remaining = [&](auto aa) -> bool {
      const int size = aa->get_npart_local();
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
    H5Part h5part("traj_reflection_dmplex_example.h5part", A, Sym<REAL>("P"),
    Sym<REAL>("V"));
    for (int stepx = 0; stepx < nsteps; stepx++) {
      nprint("step:", stepx);
      lambda_apply_timestep(static_particle_sub_group(A));
      A->hybrid_move();
      A->cell_move();

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
