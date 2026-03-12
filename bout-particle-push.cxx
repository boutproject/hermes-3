#include "bout/bout.hxx"
#include "bout/bout_types.hxx"
#include "bout/field2d.hxx"
#include "bout/output.hxx"
#include "bout/petsclib.hxx"
#include <bout/field_factory.hxx>
#include <algorithm>
#include <cmath>
#include <fmt/core.h>
#include <iostream>
#include <memory>
#include <neso_particles.hpp>
#include <neso_particles/compute_target.hpp>
#include <neso_particles/containers/cell_data.hpp>
#include <neso_particles/external_interfaces/petsc/petsc_interface.hpp>
#include <neso_particles/typedefs.hpp>
#include <netcdf>
#include <petscsystypes.h>
#include <petscviewerhdf5.h>
#include <reactions_lib/common_transformations.hpp>
#include <reactions_lib/transformation_wrapper.hpp>
#include <string>
#include <vector>
// for reactions integration
#include <reactions/reactions.hpp>
#include "include/vantage_dmplex.hxx"

#ifndef NESO_PARTICLES_PETSC
static_assert(false, "NESO-Particles was installed without PETSc support.");
#else

using namespace NESO::Particles;
using namespace VANTAGE::Reactions;

template <typename T, typename U>
inline void ASSERT_EQ(T t, U u) {
  NESOASSERT(t == u, "A check failed.");
}

void extract_ionised_density_in_place(
    Field2D& density, std::shared_ptr<CellwiseAccumulator<double>>& accumulator_transform,
    std::shared_ptr<ParticleGroup>& A_particle_group,
    std::shared_ptr<TransformationStrategy>& source_zeroer,
    std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh) {
  Mesh* bout_mesh = density.getMesh();
  std::vector<CellData<double>> accumulated_1d =
      accumulator_transform->get_cell_data("ION_SOURCE_DENSITY");
  PetscInt ic = 0;
  for (PetscInt ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (PetscInt iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      density(ix, iy) +=
          accumulated_1d[ic]->at(0, 0) / neso_mesh->dmh->get_cell_volume(ic);
      ic++;
    }
  }
  // this fills internal guards
  bout_mesh->communicate(density);
  // apply boundary conditions to fill external guards
  // density.applyBoundary();
  // extrapolate -> Neumann

  // Now set the property to zero in the accumulator buffer, ready for the next timestep.
  // Note that this does not zero the data in the original particle group
  accumulator_transform->zero_buffer("ION_SOURCE_DENSITY");
  // set the sources to zero on the particle group read for the next timestep
  source_zeroer->transform(std::make_shared<ParticleSubGroup>(A_particle_group));
}

void calculate_density_in_place(
    Field2D& density, std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG>& dg0,
    std::shared_ptr<ParticleGroup>& A_particle_group, std::vector<double>& h_project1) {
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

void calculate_fluid_moments_in_place(
    Field2D& neutral_density, Field2D& ion_density,
    std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG>& dg0,
    std::shared_ptr<ParticleGroup>& A_particle_group, std::vector<double>& h_project1,
    std::shared_ptr<CellwiseAccumulator<double>>& accumulator_transform,
    std::shared_ptr<TransformationStrategy>& source_zeroer,
    std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh) {
  // calculate all of the fluid moments required for the moment, from data in
  // the particle group. Set to zero any particle properties needed in preparation
  // for the next time step.
  calculate_density_in_place(neutral_density, dg0, A_particle_group, h_project1);
  extract_ionised_density_in_place(ion_density, accumulator_transform, A_particle_group,
                                   source_zeroer, neso_mesh);
}

double calculate_total_mass(Field2D& density,
                            std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh) {
  double local_mass = 0.0;
  double total_mass = 0.0;
  Mesh* bout_mesh = density.getMesh();
  PetscInt ic = 0;
  for (PetscInt ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (PetscInt iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      local_mass += density(ix, iy) * neso_mesh->dmh->get_cell_volume(ic);
      ic++;
    }
  }
  MPICHK(
      MPI_Allreduce(&local_mass, &total_mass, 1, MPI_DOUBLE, MPI_SUM, BoutComm::get()));
  return total_mass;
}

Options
initialise_diagnostics(Field2D& neutral_density, Field2D& ion_density,
                       std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG>& dg0,
                       std::shared_ptr<ParticleGroup>& A_particle_group,
                       std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
                       std::vector<double>& h_project1,
                       std::string particle_data_filename) {
  // Options object to use to write out diagnostic data of fluid quantities
  Options bout_output_data;
  bout_output_data["neutral_density"] = neutral_density;
  // Set the time attribute
  bout_output_data["neutral_density"].attributes["time_dimension"] = "t";
  bout_output_data["ion_density"] = ion_density;
  bout_output_data["ion_density"].attributes["time_dimension"] = "t";
  bout_output_data["total_neutral_mass"] =
      calculate_total_mass(neutral_density, neso_mesh);
  bout_output_data["total_neutral_mass"].attributes["time_dimension"] = "t";
  bout_output_data["total_ion_mass"] = calculate_total_mass(ion_density, neso_mesh);
  bout_output_data["total_ion_mass"].attributes["time_dimension"] = "t";
  Field2D total_density = ion_density + neutral_density;
  bout_output_data["total_mass"] = calculate_total_mass(total_density, neso_mesh);
  bout_output_data["total_mass"].attributes["time_dimension"] = "t";
  // std::string particle_data_filename =
  // fmt::format("bout_particle_moments_{}.nc",mpi_rank);
  bout::OptionsIO::create(particle_data_filename)->write(bout_output_data);
  return bout_output_data;
}

void update_diagnostics(Field2D& neutral_density, Field2D& ion_density,
                        std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG>& dg0,
                        std::shared_ptr<ParticleGroup>& A_particle_group,
                        std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
                        std::vector<double>& h_project1, Options& bout_output_data,
                        std::string particle_data_filename) {
  // update density in Options object and write
  bout_output_data["neutral_density"] = neutral_density;
  bout_output_data["ion_density"] = ion_density;
  Field2D total_density = ion_density + neutral_density;
  bout_output_data["total_mass"] = calculate_total_mass(total_density, neso_mesh);
  bout_output_data["total_neutral_mass"] =
      calculate_total_mass(neutral_density, neso_mesh);
  bout_output_data["total_ion_mass"] = calculate_total_mass(ion_density, neso_mesh);
  // Append data to file
  bout::OptionsIO::create({{"file", particle_data_filename}, {"append", true}})
      ->write(bout_output_data);
}

void set_initial_particle_weights(
    Field2D& initial_neutral_density,
    std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG>& dg0,
    std::shared_ptr<ParticleGroup>& A_particle_group,
    std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
    std::vector<double>& h_project1) {
  Mesh* bout_mesh = initial_neutral_density.getMesh();
  PetscInt ixy = 0;
  for (PetscInt ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (PetscInt iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      // particle_weights are copied to all particles in this cell so
      // we multiply the initial density by the volume to get particle number,
      // then divide by the number of marker particles per cell
      const REAL cell_volume = neso_mesh->dmh->get_cell_volume(ixy);
      const INT nmarkers_per_cell = A_particle_group->get_npart_cell(ixy);
      const REAL particle_weights =
          initial_neutral_density(ix, iy) * cell_volume / nmarkers_per_cell;
      h_project1.at(ixy) = particle_weights;
      ixy++;
    }
  }
  // now copy the data to internal variables
  dg0->set_dofs(1, h_project1);
  // set the data from internal variables into the weights
  dg0->evaluate(A_particle_group, Sym<REAL>("WEIGHT"));
}

void check_cell_volumes(std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
                        Mesh*& bout_mesh) {
  Coordinates* coord = bout_mesh->getCoordinates();
  PetscInt ixy = 0;
  const REAL tolerance = 1.0e-12;
  for (PetscInt ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (PetscInt iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      const REAL bout_cell_volume =
          coord->J(ix, iy) * coord->dx(ix, iy) * coord->dy(ix, iy);
      const REAL neso_cell_volume = neso_mesh->dmh->get_cell_volume(ixy);
      const bool volumes_match = (abs(bout_cell_volume - neso_cell_volume) < tolerance);
      // exit if we fail to find a match
      NESOASSERT(volumes_match,
                 fmt::format("BOUT++ mesh volume {} does not match NESO-Particles mesh "
                             "volume {} for ix = {} iy = {} \n Ignore this message by "
                             "setting [neso_particles] test_cell_volumes = false",
                             bout_cell_volume, neso_cell_volume, ix, iy));
      ixy++;
    }
  }
}

REAL cell_length(std::vector<std::vector<REAL>> cell_vertices, INT iv1, INT iv2, INT iv3, INT iv4){
  const REAL Rlength2 = std::pow(0.5*(cell_vertices.at(iv1).at(0) + cell_vertices.at(iv2).at(0) -
                 cell_vertices.at(iv3).at(0) - cell_vertices.at(iv4).at(0)),2.0);
  const REAL Zlength2 = std::pow(0.5*(cell_vertices.at(iv1).at(1) + cell_vertices.at(iv2).at(1) -
                 cell_vertices.at(iv3).at(1) - cell_vertices.at(iv4).at(1)),2.0);
  REAL length = std::pow(Zlength2+Rlength2,0.5);
  return length;
}


void check_cell_centres(std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
                        Mesh*& bout_mesh, BoutReal absolute_tolerance, BoutReal relative_tolerance) {
  Coordinates* coord = bout_mesh->getCoordinates();
  // get (R,Z) of cell centres in Hypnotoad grid
  Field2D Rxy;
  Field2D Zxy;
  bout_mesh->get(Rxy, "Rxy");
  bout_mesh->get(Zxy, "Zxy");
  // compare to cell centres calculated from cell corners
  std::vector<std::vector<REAL>> cell_vertices;
  PetscInt ixy = 0;
  for (PetscInt ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (PetscInt iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      const REAL bout_Rxy = Rxy(ix, iy);
      const REAL bout_Zxy = Zxy(ix, iy);
      neso_mesh->dmh->get_cell_vertices(ixy, cell_vertices);

      REAL neso_Rxy = 0.0;
      REAL neso_Zxy = 0.0;
      for (PetscInt iv = 0; iv < 4; iv++){
        neso_Rxy += cell_vertices.at(iv).at(0);
        neso_Zxy += cell_vertices.at(iv).at(1);
      }
      neso_Rxy /= 4.0;
      neso_Zxy /= 4.0;
      // get lengths of cell across the two dimensions
      const REAL cell_length_a = cell_length(cell_vertices, 0, 1, 2, 3);
      const REAL cell_length_b = cell_length(cell_vertices, 0, 3, 2, 1);
      const REAL min_cell_length = std::min(cell_length_a, cell_length_b);
      // we compare the difference in cell centres to the absolute tolerance and
      // the relative tolerance formed by comparing to the smallest length across the cell
      const REAL tolerance = absolute_tolerance + min_cell_length*relative_tolerance;
      const bool centres_match = (abs(neso_Rxy - bout_Rxy) < tolerance) &&
                                 (abs(neso_Zxy - bout_Zxy) < tolerance);
      // exit if we fail to find a match
      NESOASSERT(centres_match,
                 fmt::format("Hypnotoad/BOUT++ cell centre (R, Z) ({}, {}) does not match NESO-Particles mesh "
                             "cell centre ({}, {}) for ix = {} iy = {} \n"
                             "The cell height and width are {} {} \n"
                             "The displacements in R and Z are {} {} \n"
                             "Ignore this message by "
                             "setting [neso_particles] test_cell_centres = false\n Relax the "
                             "tolerance used in this check by increasing\n"
                             "[neso_particles] cell_centre_absolute_tolerance = {}\n"
                             "[neso_particles] cell_centre_relative_tolerance = {}",
                             bout_Rxy, bout_Zxy, neso_Rxy, neso_Zxy, ix, iy,
                             cell_length_a, cell_length_b, abs(neso_Rxy - bout_Rxy), abs(neso_Zxy - bout_Zxy),
                             absolute_tolerance, relative_tolerance));
      ixy++;
    }
  }
}

void check_mass_conservation(double total_mass_final, double total_mass_initial,
                double remove_threshold) {
  double rtol = 1.0e-13;
  double atol = 1.0e-13;
  bool mass_conserved = (abs(total_mass_final - total_mass_initial) < rtol*total_mass_initial + atol);
  // exit if we fail to find conservation
  NESOASSERT(mass_conserved,
              fmt::format("Initial total mass {} does not match "
                          "final total mass {} \n Ignore this message by "
                          "setting [neso_particles] test_mass_conservation = false",
                          total_mass_initial, total_mass_final));
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
  Field2D initial_neutral_density{bout_mesh};
  initial_neutral_density = opt["mesh"]["initial_neutral_density"].as<Field2D>();
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
    Field2D ion_density = Field2D(0.0, bout_mesh);
    Field2D neutral_density = Field2D(0.0, bout_mesh);
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
    if (Options::root()["neso_particles"]["test_cell_volumes"].withDefault(true)) {
      check_cell_volumes(neso_mesh, bout_mesh);
    }
    if (Options::root()["neso_particles"]["test_cell_centres"].withDefault(true)){
      check_cell_centres(neso_mesh,bout_mesh,
        Options::root()["neso_particles"]["cell_centre_absolute_tolerance"].withDefault(1.0e-12),
        Options::root()["neso_particles"]["cell_centre_relative_tolerance"].withDefault(0.0));
    }
    // create a Reactions particle spec
    auto particle_spec_builder = ParticleSpecBuilder(ndim);
    auto electron_species = Species("ELECTRON");
    auto main_species = Species("ION", 1.0, 0.0, 0);
    std::vector<Species> fluid_species = {electron_species, main_species};
    particle_spec_builder.add_particle_prop(Properties<REAL>(
        fluid_species,
        std::vector<int>{default_properties.temperature, default_properties.density,
                         default_properties.source_energy,
                         default_properties.source_density}));
    particle_spec_builder.add_particle_prop(
        Properties<REAL>(fluid_species,
                         std::vector<int>{default_properties.source_momentum}),
        ndim);
    ParticleSpec additional_props{ParticleProp(Sym<REAL>("TSP"), 2)};
    particle_spec_builder.add_particle_spec(additional_props);
    ParticleSpec particle_spec = particle_spec_builder.get_particle_spec();

    // Create a Particle group with our specied particle properties.
    auto A_particle_group =
        std::make_shared<ParticleGroup>(domain, particle_spec, sycl_target);

    // Create some particle data
    const int mpi_rank = sycl_target->comm_pair.rank_parent;
    std::mt19937 rng_pos(52234234 + mpi_rank);
    std::mt19937 rng_vel(52234231 + mpi_rank);
    std::vector<std::vector<double>> positions;
    std::vector<int> particle_cell_ids;

    uniform_within_dmplex_cells(neso_mesh, npart_per_cell, positions, particle_cell_ids,
                                &rng_pos);

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
      initial_distribution[Sym<REAL>("ION_SOURCE_DENSITY")][px][0] = 0.0;
      initial_distribution[Sym<REAL>("ION_SOURCE_ENERGY")][px][0] = 0.0;
      initial_distribution[Sym<REAL>("ELECTRON_DENSITY")][px][0] = 1.0;
      initial_distribution[Sym<REAL>("ELECTRON_SOURCE_DENSITY")][px][0] = 0.0;
      initial_distribution[Sym<REAL>("ELECTRON_SOURCE_ENERGY")][px][0] = 0.0;
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
    const REAL remove_threshold =
        Options::root()["VANTAGE_reactions"]["remove_threshold"].withDefault(1.0e-10);
    const REAL merge_threshold =
        Options::root()["VANTAGE_reactions"]["merge_threshold"].withDefault(1.0e-2);

    auto remove_wrapper = std::make_shared<TransformationWrapper>(
        std::vector<std::shared_ptr<MarkingStrategy>>{
            make_marking_strategy<ComparisonMarkerSingle<REAL, LessThanComp>>(
                Sym<REAL>("WEIGHT"), remove_threshold)},
        std::dynamic_pointer_cast<TransformationStrategy>(
            std::make_shared<SimpleRemovalTransformationStrategy>()));

    auto merge_wrapper = std::make_shared<TransformationWrapper>(
        std::vector<std::shared_ptr<MarkingStrategy>>{
            make_marking_strategy<ComparisonMarkerSingle<REAL, LessThanComp>>(
                Sym<REAL>("WEIGHT"), merge_threshold)},
        make_transformation_strategy<MergeTransformationStrategy<ndim>>());

    auto accumulator_transform = std::make_shared<CellwiseAccumulator<REAL>>(
        A_particle_group, std::vector<std::string>{"ION_SOURCE_DENSITY"});
    // std::vector<CellData<double>> accumulated_1d =
    // accumulator_transform->get_cell_data("ION_SOURCE_DENSITY");
    auto accumulator_real_transform_wrapper = std::make_shared<TransformationWrapper>(
        std::dynamic_pointer_cast<TransformationStrategy>(accumulator_transform));
    auto source_zeroer = make_transformation_strategy<ParticleDatZeroer<REAL>>(
        std::vector<std::string>{"ION_SOURCE_DENSITY"});

    std::vector<std::shared_ptr<TransformationWrapper>> child_transforms =
        std::vector{merge_wrapper, remove_wrapper};
    std::vector<std::shared_ptr<TransformationWrapper>> parent_transforms =
        std::vector{accumulator_real_transform_wrapper, merge_wrapper, remove_wrapper};

    auto reaction_controller = ReactionController(parent_transforms, child_transforms);
    // add ionisation to the controller
    reaction_controller.add_reaction(
        std::make_shared<decltype(ionisation_reaction)>(ionisation_reaction));

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
        reflection->execute(gx.second, Sym<REAL>("POSITION"), Sym<REAL>("VELOCITY"),
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
          [=](auto VELOCITY, auto POSITION, auto TSP) {
            const REAL dt_left = dt - TSP.at(0);
            if (dt_left > 0.0) {
              POSITION.at(0) += dt_left * VELOCITY.at(0);
              POSITION.at(1) += dt_left * VELOCITY.at(1);
              TSP.at(0) = dt;
              TSP.at(1) = dt_left;
            }
          },
          Access::read(Sym<REAL>("VELOCITY")), Access::write(Sym<REAL>("POSITION")),
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
    H5Part h5part("traj_reflection_dmplex_example.h5part", A_particle_group,
                  Sym<REAL>("POSITION"), Sym<REAL>("VELOCITY"));

    // allocate buffer vector for scalar projection/evaluation of NESO-Particles
    // properties
    std::vector<REAL> h_project1(num_cells_owned);
    // set weights from a Field2D from BOUT
    set_initial_particle_weights(initial_neutral_density, dg0, A_particle_group,
                                 neso_mesh, h_project1);
    // update fluid moments
    calculate_fluid_moments_in_place(neutral_density, ion_density, dg0, A_particle_group,
                                     h_project1, accumulator_transform, source_zeroer,
                                     neso_mesh);
    // diagnose the initial condition
    std::string particle_data_filename =
        fmt::format("bout_particle_moments_{}.nc", mpi_rank);
    Options bout_output_data =
        initialise_diagnostics(neutral_density, ion_density, dg0, A_particle_group,
                               neso_mesh, h_project1, particle_data_filename);
    // mass for conservation check
    Field2D total_density = neutral_density + ion_density;
    double total_mass_initial = calculate_total_mass(total_density, neso_mesh);
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
      // A_particle_group->print(Sym<REAL>("POSITION"), Sym<INT>("ID"),
      // Sym<REAL>("WEIGHT"), Sym<REAL>("ION_SOURCE_DENSITY")); update fluid moments
      calculate_fluid_moments_in_place(neutral_density, ion_density, dg0,
                                       A_particle_group, h_project1,
                                       accumulator_transform, source_zeroer, neso_mesh);
      // diagnose timestep stepx
      update_diagnostics(neutral_density, ion_density, dg0, A_particle_group, neso_mesh,
                         h_project1, bout_output_data, particle_data_filename);
    }
    // uncomment to write a trajectory
    h5part.close();

    // mass for conservation check
    total_density = neutral_density + ion_density;
    double total_mass_final = calculate_total_mass(total_density, neso_mesh);
    if (Options::root()["neso_particles"]["test_mass_conservation"].withDefault(true)) {
      check_mass_conservation(total_mass_final, total_mass_initial, remove_threshold);
    }
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
