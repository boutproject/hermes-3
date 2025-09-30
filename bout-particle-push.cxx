#include "bout/petsclib.hxx"
#include <neso_particles.hpp>
#include <neso_particles/external_interfaces/petsc/petsc_interface.hpp>
#include "bout/bout.hxx"

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
  PETSCCHK(PetscInitializeNoArguments());
  auto sycl_target = std::make_shared<SYCLTarget>(0, PETSC_COMM_WORLD);
  const int mpi_size = sycl_target->comm_pair.size_parent;
  const int mpi_rank = sycl_target->comm_pair.rank_parent;

  // First we setup the topology of the mesh.
  PetscInt num_cells_owned = mpi_size;

  std::vector<PetscInt> cells;
  cells.reserve(num_cells_owned * 4);

  // We are careful to list the vertices in counter clock-wise order, this might
  // matter.
  for (PetscInt cx = 0; cx < num_cells_owned; cx++) {
    // These are global indices not local indices.
    PetscInt vertex_sw = cx + mpi_rank * (mpi_size + 1);
    PetscInt vertex_se = vertex_sw + 1;
    PetscInt vertex_ne = vertex_se + mpi_size + 1;
    PetscInt vertex_nw = vertex_ne - 1;
    cells.push_back(vertex_sw);
    cells.push_back(vertex_se);
    cells.push_back(vertex_ne);
    cells.push_back(vertex_nw);
  }

  /*
   * Each rank owns a contiguous block of global indices. We label our indices
   * lexicographically (row-wise). Sorting out the global vertex indexing is
   * probably one of the more tedious parts.
   */
  PetscInt num_vertices_owned = mpi_size + 1;
  if (mpi_rank == mpi_size - 1) {
    // Make the last rank own the top edge.
    num_vertices_owned += mpi_size + 1;
  }

  /*
   * Create the coordinates for the block of vertices we pass to petsc. For an
   * existing mesh in memory this step will probably involve some MPI
   * communication to gather the blocks of coordinates on the ranks which pass
   * them to PETSc.
   */
  std::vector<PetscScalar> vertex_coords(num_vertices_owned * 2);
  for (int px = 0; px < mpi_size + 1; px++) {
    // our cell extent is 1.0
    vertex_coords.at(px * 2 + 0) = px;
    vertex_coords.at(px * 2 + 1) = mpi_rank;
  }
  if (mpi_rank == mpi_size - 1) {
    for (int px = 0; px < mpi_size + 1; px++) {
      // our cell extent is 1.0
      vertex_coords.at((mpi_size + 1 + px) * 2 + 0) = px;
      vertex_coords.at((mpi_size + 1 + px) * 2 + 1) = mpi_rank + 1;
    }
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

  // Label subsections of the boundary by specifing pairs of vertices and using
  // the label_dmplex_edges helper function.
  std::vector<PetscInt> vertex_starts, vertex_ends, edge_labels;

  if (mpi_rank == mpi_size - 1) {
    // Top edge
    for (int px = 0; px < mpi_size; px++) {
      const PetscInt tx = (mpi_size + 1) * mpi_size + px;
      vertex_starts.push_back(tx);
      vertex_ends.push_back(tx + 1);
      // Label the top edge with label 200
      edge_labels.push_back(200);
    }
  }

  PetscInterface::label_dmplex_edges(dm, PetscInterface::face_sets_label,
                                     vertex_starts, vertex_ends, edge_labels);

  /*
   *
   *
   *
   *
   *
   * Below here is testing of the DMPlex
   *
   *
   *
   *
   *
   *
   */

  // Create the map from local point numbering to global point numbering.
  IS global_point_numbers;
  PETSCCHK(DMPlexCreatePointNumbering(dm, &global_point_numbers));

  // Get the DMPlex point start/end for the cells this process owns in the
  // DMPlex
  PetscInt cell_local_start, cell_local_end;
  PETSCCHK(DMPlexGetDepthStratum(dm, 2, &cell_local_start, &cell_local_end));
  // Map the local cell indices to global indices
  PetscInt point_start, point_end;
  PETSCCHK(DMPlexGetChart(dm, &point_start, &point_end));
  const PetscInt *ptr;
  PETSCCHK(ISGetIndices(global_point_numbers, &ptr));
  PetscInt local_index = 0;
  PetscInt correct_cell_index_start = mpi_rank * mpi_size;
  for (PetscInt point = cell_local_start; point < cell_local_end; point++) {
    PetscInt cell_index = ptr[point - point_start];
    ASSERT_EQ(correct_cell_index_start + local_index, cell_index);
    local_index++;
  }
  PETSCCHK(ISRestoreIndices(global_point_numbers, &ptr));

  // Check the vertices of the local cells are the vertices we expect.
  // This also checks that the ordering and parallel decomposition is what we
  // expect
  PetscInterface::DMPlexHelper dmh(PETSC_COMM_WORLD, dm);
  std::vector<std::vector<REAL>> vertices;
  std::set<std::vector<REAL>> vertices_to_test;
  std::set<std::vector<REAL>> vertices_correct;

  for (int cx = 0; cx < mpi_size; cx++) {
    dmh.get_cell_vertices(cx, vertices);
    vertices_to_test.clear();
    vertices_correct.clear();

    for (auto vx : vertices) {
      vertices_to_test.insert(vx);
    }

    vertices_correct.insert(
        {static_cast<REAL>(cx), static_cast<REAL>(mpi_rank)});
    vertices_correct.insert(
        {static_cast<REAL>(cx + 1), static_cast<REAL>(mpi_rank)});
    vertices_correct.insert(
        {static_cast<REAL>(cx), static_cast<REAL>(mpi_rank + 1)});
    vertices_correct.insert(
        {static_cast<REAL>(cx + 1), static_cast<REAL>(mpi_rank + 1)});
    ASSERT_EQ(vertices_correct, vertices_to_test);
  }

  // Check the edges are labelled correctly
  {
    auto lambda_is_on_boundary = [&](auto coords) -> bool {
      auto lambda_check_coord = [&](auto coord) -> bool {
        if (std::abs(coord[0]) < 1.0e-14) {
          return true;
        }
        if (std::abs(coord[0] - (mpi_size)) < 1.0e-14) {
          return true;
        }
        if (std::abs(coord[1]) < 1.0e-14) {
          return true;
        }
        if (std::abs(coord[1] - (mpi_size)) < 1.0e-14) {
          return true;
        }
        return false;
      };
      return lambda_check_coord(coords) && lambda_check_coord(coords + 2);
    };

    PetscInt edge_start, edge_end;
    PETSCCHK(DMPlexGetDepthStratum(dm, 1, &edge_start, &edge_end));
    DMLabel label;
    PETSCCHK(DMGetLabel(dm, PetscInterface::face_sets_label, &label));
    for (PetscInt edge = edge_start; edge < edge_end; edge++) {
      PetscInt size;
      PETSCCHK(DMPlexGetConeSize(dm, edge, &size));
      ASSERT_EQ(size, 2);

      PetscBool is_dg;
      PetscInt num_coords;
      const PetscScalar *array;
      PetscScalar *coords;
      PETSCCHK(DMPlexGetCellCoordinates(dm, edge, &is_dg, &num_coords, &array,
                                        &coords));

      if (lambda_is_on_boundary(coords)) {
        PetscInt value;
        PETSCCHK(DMLabelGetValue(label, edge, &value));
        // is this a top edge, we labelled all the top edges 200
        if ((std::abs(coords[0] - mpi_size) < 1.0e-14) &&
            (std::abs(coords[1] - mpi_size) < 1.0e-14) &&
            (std::abs(coords[2] - mpi_size) < 1.0e-14) &&
            (std::abs(coords[3] - mpi_size) < 1.0e-14)) {
          ASSERT_EQ(value, 200);
        } else {
          ASSERT_EQ(value, 100);
        }
      }

      PETSCCHK(DMPlexRestoreCellCoordinates(dm, edge, &is_dg, &num_coords,
                                            &array, &coords));
    }
  }

  dmh.free();

  PETSCCHK(ISDestroy(&global_point_numbers));

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
    auto mesh = std::make_shared<PetscInterface::DMPlexInterface>(
        dm, 0, MPI_COMM_WORLD);
    // Create a mapper for mapping particles into cells.
    auto mapper =
        std::make_shared<PetscInterface::DMPlexLocalMapper>(sycl_target, mesh);
    // Create a domain from the mesh and the mapper.
    auto domain = std::make_shared<Domain>(mesh, mapper);

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

    uniform_within_dmplex_cells(mesh, npart_per_cell, positions, cells,
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
        sycl_target, mesh, boundary_groups);
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
    // H5Part h5part("traj_reflection_dmplex_example.h5part", A, Sym<REAL>("P"),
    // Sym<REAL>("V"));
    for (int stepx = 0; stepx < nsteps; stepx++) {
      nprint("step:", stepx);
      lambda_apply_timestep(static_particle_sub_group(A));
      A->hybrid_move();
      A->cell_move();

      // uncomment to write a trajectory
      // h5part.write();
    }

    // uncomment to write a trajectory
    // h5part.close();

    // Boundary interaction objects require a free call.
    b2d->free();
    // NESO-Particles mesh objects must have free called on them.
    mesh->free();
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
