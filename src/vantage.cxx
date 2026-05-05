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
#include "../include/vantage.hxx"
#include "../include/vantage_dmplex.hxx"
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

// Make path for any output file
std::string make_output_path(const std::string& filename, Options& alloptions) {
  std::string output_dir = alloptions["restart_files"]["path"].withDefault(
      alloptions["datadir"].withDefault<std::string>("data"));
  return fmt::format("{}/{}", output_dir, filename);
}

void calculate_neutral_density_in_place(
    Field2D& density, std::shared_ptr<PetscInterface::DMPlexProjectEvaluateDG>& dg0,
    std::shared_ptr<ParticleGroup>& A_particle_group, std::vector<double>& h_project1) {
  Mesh* bout_mesh = density.getMesh();
  // get a density by projecting the particle property WEIGHT to the bout_mesh
  dg0->project(A_particle_group, Sym<REAL>("WEIGHT"));
  // std::vector<REAL> h_project1;
  dg0->get_dofs(1, h_project1);
  std::size_t ic = 0;
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
initialise_diagnostics(Mesh* bout_mesh, Field2D& neutral_density, Field2D& ion_density,
                       std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
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
  bout_output_data["t_array"] = 0.0;
  bout_output_data["t_array"].attributes["time_dimension"] = "t";
  // Mesh metadata
  bout_mesh->outputVars(bout_output_data);

  bout::OptionsIO::create(particle_data_filename)->write(bout_output_data);
  return bout_output_data;
}

void update_diagnostics(Field2D& neutral_density, Field2D& ion_density,
                        std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
                        Options& bout_output_data, std::string particle_data_filename,
                        BoutReal particle_time) {
  // update density in Options object and write
  bout_output_data["neutral_density"] = neutral_density;
  bout_output_data["ion_density"] = ion_density;
  Field2D total_density = ion_density + neutral_density;
  bout_output_data["total_mass"] = calculate_total_mass(total_density, neso_mesh);
  bout_output_data["total_neutral_mass"] =
      calculate_total_mass(neutral_density, neso_mesh);
  bout_output_data["total_ion_mass"] = calculate_total_mass(ion_density, neso_mesh);
  bout_output_data["t_array"] = particle_time;
  // bout_output_data["t_array"] = 0.0;
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
      const REAL cell_volume = neso_mesh->dmh->get_cell_volume(static_cast<int>(ixy));
      const INT nmarkers_per_cell =
          A_particle_group->get_npart_cell(static_cast<int>(ixy));
      const REAL particle_weights = initial_neutral_density(ix, iy) * cell_volume
                                    / static_cast<BoutReal>(nmarkers_per_cell);
      h_project1.at(static_cast<std::size_t>(ixy)) = particle_weights;
      ixy++;
    }
  }
  // now copy the data to internal variables
  dg0->set_dofs(1, h_project1);
  // set the data from internal variables into the weights
  dg0->evaluate(A_particle_group, Sym<REAL>("WEIGHT"));
}

void check_cell_volumes(std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
                        Mesh*& bout_mesh, Options& alloptions) {
  Coordinates* coord = bout_mesh->getCoordinates();
  PetscInt ixy = 0;
  const REAL tolerance = 1.0e-12;
  for (PetscInt ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (PetscInt iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {

      // Convert to SI: dx is m^2 T, J is m/T, dy is unitless
      // so J * dx * dy = m^3
      const BoutReal meters = get<BoutReal>(alloptions["units"]["meters"]);
      const BoutReal meters_cubed = meters * meters * meters;

      const BoutReal bout_cell_volume =
          coord->J(ix, iy) * coord->dx(ix, iy) * coord->dy(ix, iy) * meters_cubed;

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

REAL cell_length(std::vector<std::vector<REAL>>& cell_vertices, std::size_t iv1,
                 std::size_t iv2, std::size_t iv3, std::size_t iv4) {
  const REAL Rlength2 =
      std::pow(0.5
                   * (cell_vertices.at(iv1).at(0) + cell_vertices.at(iv2).at(0)
                      - cell_vertices.at(iv3).at(0) - cell_vertices.at(iv4).at(0)),
               2.0);
  const REAL Zlength2 =
      std::pow(0.5
                   * (cell_vertices.at(iv1).at(1) + cell_vertices.at(iv2).at(1)
                      - cell_vertices.at(iv3).at(1) - cell_vertices.at(iv4).at(1)),
               2.0);
  REAL length = std::pow(Zlength2 + Rlength2, 0.5);
  return length;
}

void check_cell_centres(std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh,
                        Mesh*& bout_mesh, BoutReal absolute_tolerance,
                        BoutReal relative_tolerance) {
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
      for (std::size_t iv = 0; iv < 4; iv++) {
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
      const REAL tolerance = absolute_tolerance + min_cell_length * relative_tolerance;
      const bool centres_match = (abs(neso_Rxy - bout_Rxy) < tolerance)
                                 && (abs(neso_Zxy - bout_Zxy) < tolerance);
      // exit if we fail to find a match
      NESOASSERT(
          centres_match,
          fmt::format("Hypnotoad/BOUT++ cell centre (R, Z) ({}, {}) does not match "
                      "NESO-Particles mesh "
                      "cell centre ({}, {}) for ix = {} iy = {} \n"
                      "The cell height and width are {} {} \n"
                      "The displacements in R and Z are {} {} \n"
                      "Ignore this message by "
                      "setting [neso_particles] test_cell_centres = false\n Relax the "
                      "tolerance used in this check by increasing\n"
                      "[neso_particles] cell_centre_absolute_tolerance = {}\n"
                      "[neso_particles] cell_centre_relative_tolerance = {}",
                      bout_Rxy, bout_Zxy, neso_Rxy, neso_Zxy, ix, iy, cell_length_a,
                      cell_length_b, abs(neso_Rxy - bout_Rxy), abs(neso_Zxy - bout_Zxy),
                      absolute_tolerance, relative_tolerance));
      ixy++;
    }
  }
}

void check_mass_conservation(double total_mass_final, double total_mass_initial) {
  double rtol = 1.0e-13;
  double atol = 1.0e-13;
  bool mass_conserved =
      (abs(total_mass_final - total_mass_initial) < rtol * total_mass_initial + atol);
  // exit if we fail to find conservation
  NESOASSERT(mass_conserved,
             fmt::format("Initial total mass {} does not match "
                         "final total mass {} \n Ignore this message by "
                         "setting [neso_particles] test_mass_conservation = false",
                         total_mass_initial, total_mass_final));
}

VantageSourceManager::VantageSourceManager(
    std::shared_ptr<PetscInterface::DMPlexInterface>& neso_mesh, Mesh* bout_mesh,
    Options& units)
    : bout_mesh(bout_mesh), neso_mesh(neso_mesh), units(units) {}

// Register new source with the manager and initialise its data
void VantageSourceManager::add_source(
    const std::string& hermes_source_name, const std::string& vantage_source_name,
    std::shared_ptr<CellwiseAccumulator<REAL>> accumulator,
    std::shared_ptr<ParticleGroup> particle_group,
    std::shared_ptr<TransformationStrategy> zeroer) {

  Field2D source_data{bout_mesh};
  source_data = 0.0;

  VantageSource source{
      hermes_source_name, vantage_source_name, accumulator, particle_group, zeroer,
      source_data};

  this->sources[hermes_source_name] = source;
}

// Return source data
Field2D VantageSourceManager::get_data(const std::string& hermes_source_name) {
  return this->sources[hermes_source_name].source_data;
}

// Update the source from VANTAGE and reset the VANTAGE data/accumulator
void VantageSourceManager::update_source(const std::string& hermes_source_name,
                                         double dt) {
  VantageSource& source = this->sources[hermes_source_name];

  std::vector<CellData<double>> accumulated_1d =
      source.accumulator->get_cell_data(source.vantage_source_name);

  std::size_t ic = 0;
  for (int ix = bout_mesh->xstart; ix <= bout_mesh->xend; ix++) {
    for (int iy = bout_mesh->ystart; iy <= bout_mesh->yend; iy++) {
      source.source_data(ix, iy) =
          accumulated_1d[ic]->at(0, 0)
          / neso_mesh->dmh->get_cell_volume(static_cast<int>(ic));
      ic++;
    }
  }
  // Divite by dt so that the source is in per second
  // TODO: Ensure this has correct units and units consistent with neutral density
  source.source_data /= dt;

  // Fill internal guards
  bout_mesh->communicate(source.source_data);
  // Reset the accumulator object
  source.accumulator->zero_buffer(source.vantage_source_name);
  // Reset the accumulated source data on the particle
  source.zeroer->transform(std::make_shared<ParticleSubGroup>(source.particle_group));
}

// Update all sources
void VantageSourceManager::update_all_sources(double dt) {
  for (auto& [hermes_source_name, source] : this->sources) {
    update_source(hermes_source_name, dt);
  }
}

Vantage::Vantage(std::string name, Options& alloptions, Solver* UNUSED(solver))
    : Component({readOnly("species:d+:density", Regions::Interior),
                 readWrite("species:d+:density")}) {

  // TODO: Put proper permissions in

  // int main(int argc, char** argv) {
  // initialise_mpi(&argc, &argv);
  // attempt to call BOUT to
  // get information from a BOUT
  // mesh object
  // N.B. Comment the next three lines
  // to permit compilation as is

  // BoutInitialise(argc, argv);
  // Mesh* bout_mesh = Mesh::create(&Options::root()["mesh"]);
  // TODO: tidy up the above

  Options& options = alloptions[name];        // [vantage]
  Options& mesh_options = alloptions["dmplex"]; // [mesh]
  Options& units = alloptions["units"];

  BoutReal N_w =
      units["N_w"]
          .doc("Normalisation parameter: neutral particle density per unit weight")
          .withDefault(2.0);
  BoutReal Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
  BoutReal Tnorm = get<BoutReal>(units["eV"]);

  Mesh* bout_mesh = bout::globals::mesh;
  sycl_target = std::make_shared<SYCLTarget>(0, BoutComm::get());
  // keep dmplex_h5_filename in vantage.cxx to retain access to make_output_path()
  // which should presumably not need to exist within the hermes-3 library
  std::string dmplex_h5_filename = mesh_options["dmplex_h5_filename"]
                                       .doc("Filename to use for saving the DMPlex mesh")
                                       .withDefault("hypnotoad_dmplex_mesh_output.h5");
  dm = create_dmplex_from_Bout_mesh(bout_mesh, mesh_options, sycl_target,
                                    make_output_path(dmplex_h5_filename, alloptions));
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

    // Normalisations
    //  TODO: Implement real normalisation consistent with Hermes-3
    //  TODO: Normalise DMPlex

    // Initial neutral parameters
    Field2D initial_neutral_density{bout_mesh};
    initial_neutral_density =
        options["initial_neutral_density"]
            .doc(
                "Initial neutral density for VANTAGE kinetic neutrals (normalised units)")
            .as<Field2D>();
    const int npart_per_cell = options["npart_per_cell"]
                                   .doc("Number of VANTAGE kinetic neutral particles per "
                                        "cell during initialisation")
                                   .withDefault(1);

    // Plasma parameters
    const BoutReal background_ion_temperature =
        options["background_ion_temperature"].withDefault(1.0);
    const BoutReal background_ion_density =
        options["background_ion_density"].withDefault(1.0);
    const BoutReal background_ion_Vx = options["background_ion_Vx"].withDefault(0.0);
    const BoutReal background_ion_Vy = options["background_ion_Vy"].withDefault(0.0);
    const std::vector<BoutReal> V_background = {background_ion_Vx, background_ion_Vy};

    // Reaction settings
    const REAL iz_rate = options["iz_rate"].withDefault(1.0);
    const REAL rec_rate = options["rec_rate"].withDefault(1.0);
    const int rec_markers_per_cell = options["rec_markers_per_cell"].withDefault(1000);

    // Other settings
    const int ndim = 2;
    const REAL dt = options["dt"]
                        .doc("Timestep to use for VANTAGE kinetic neutrals")
                        .withDefault(0.01);
    const int nsteps = options["nsteps"]
                           .doc("Number of timesteps to use for VANTAGE kinetic neutrals")
                           .withDefault(10);
    const int rng_samples = options["rng_samples"]
                                .doc("Number of RNG samples to prepare per-particle")
                                .withDefault(40);

    BoutReal particle_time = 0.0;
    ion_density = Field2D(background_ion_density, bout_mesh);
    neutral_density = Field2D(0.0, bout_mesh);
    // Create a mesh interface from the DM
    neso_mesh = std::make_shared<PetscInterface::DMPlexInterface>(dm, 0, BoutComm::get());
    // Create a mapper for mapping particles into cells.
    auto mapper =
        std::make_shared<PetscInterface::DMPlexLocalMapper>(sycl_target, neso_mesh);
    // Create a domain from the neso_mesh and the mapper.
    auto domain = std::make_shared<Domain>(neso_mesh, mapper);
    // Get the number of cells in the mesh owned on this process
    int num_cells_owned = neso_mesh->get_cell_count();
    // if requested, check that neso_mesh cell volumes are identical
    // to bout_mesh cell volumes, otherwise, exit.
    if (mesh_options["test_dmplex_cell_volumes"].withDefault(true)) {
      check_cell_volumes(neso_mesh, bout_mesh, alloptions);
    }
    if (mesh_options["test_dmplex_cell_centres"].withDefault(true)) {
      check_cell_centres(
          neso_mesh, bout_mesh,
          mesh_options["dmplex_cell_centre_absolute_tolerance"].withDefault(1.0e-12),
          mesh_options["dmplex_cell_centre_relative_tolerance"].withDefault(0.0));
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
    ParticleSpec additional_props{ParticleProp(Sym<REAL>("TSP"), 2),
                                  ParticleProp(Sym<REAL>("FLUID_DENSITY"), 1),
                                  ParticleProp(Sym<REAL>("FLUID_FLOW_SPEED"), ndim),
                                  ParticleProp(Sym<REAL>("FLUID_TEMPERATURE"), 1),
                                  ParticleProp(Sym<INT>("N_CELL"), 1)};

    particle_spec_builder.add_particle_spec(additional_props);
    ParticleSpec particle_spec = particle_spec_builder.get_particle_spec();

    // Create a Particle group with our specied particle properties.
    auto A_particle_group =
        std::make_shared<ParticleGroup>(domain, particle_spec, sycl_target);

    // Create some particle data
    const int mpi_rank = sycl_target->comm_pair.rank_parent;
    std::mt19937 rng_pos(static_cast<std::mt19937::result_type>(52234234 + mpi_rank));
    std::mt19937 rng_vel(static_cast<std::mt19937::result_type>(52234231 + mpi_rank));
    std::vector<std::vector<double>> positions;
    std::vector<int> particle_cell_ids;

    uniform_within_dmplex_cells(neso_mesh, npart_per_cell, positions, particle_cell_ids,
                                &rng_pos);

    const int N_actual = static_cast<int>(particle_cell_ids.size());
    auto velocities =
        NESO::Particles::normal_distribution(N_actual, 2, 0.0, 1.0, rng_vel);

    int id_offset = 0;
    MPICHK(MPI_Exscan(&N_actual, &id_offset, 1, MPI_INT, MPI_SUM,
                      sycl_target->comm_pair.comm_parent));

    // This is host space to create particle data in before pushing the
    // particles into the ParticleGroup
    ParticleSet initial_distribution(N_actual, particle_spec);
    for (int px = 0; px < N_actual; px++) {
      const auto pxu = static_cast<std::size_t>(px);

      for (int dimx = 0; dimx < ndim; dimx++) {
        const auto dimu = static_cast<std::size_t>(dimx);

        initial_distribution[Sym<REAL>("POSITION")][px][dimx] = positions[dimu][pxu];
        initial_distribution[Sym<REAL>("VELOCITY")][px][dimx] = velocities[dimu][pxu];
      }
      initial_distribution[Sym<REAL>("ION_DENSITY")][px][0] = background_ion_density;
      initial_distribution[Sym<REAL>("ION_SOURCE_DENSITY")][px][0] = 0.0;
      initial_distribution[Sym<REAL>("ION_SOURCE_ENERGY")][px][0] = 0.0;
      initial_distribution[Sym<REAL>("ELECTRON_DENSITY")][px][0] = background_ion_density;
      initial_distribution[Sym<REAL>("ELECTRON_SOURCE_DENSITY")][px][0] = 0.0;
      initial_distribution[Sym<REAL>("ELECTRON_SOURCE_ENERGY")][px][0] = 0.0;
      for (int dimx = 0; dimx < ndim; dimx++) {
        initial_distribution[Sym<REAL>("ION_SOURCE_MOMENTUM")][px][dimx] = 0.0;
        initial_distribution[Sym<REAL>("ELECTRON_SOURCE_MOMENTUM")][px][dimx] = 0.0;
      }
      initial_distribution[Sym<INT>("CELL_ID")][px][0] = particle_cell_ids.at(pxu);
      initial_distribution[Sym<INT>("ID")][px][0] = px + id_offset;
      initial_distribution[Sym<REAL>("WEIGHT")][px][0] = 1.0;
    }
    // Add the new particles to the particle group
    A_particle_group->add_particles_local(initial_distribution);
    // make pointer to projection object
    auto dg0 = std::make_shared<PetscInterface::DMPlexProjectEvaluateDG>(
        neso_mesh, sycl_target, "DG", 0);

    // RNG kernel
    // Used for sampling from velocity distribution for REC/CX
    // ------------------------------------------------------------------------------

    auto rng_kernel =
        get_uniform_rng_kernel(sycl_target, static_cast<size_t>(rng_samples));

    // Ionisation reaction
    // ------------------------------------------------------------------------------

    auto iz_rate_data = FixedRateData(iz_rate);
    main_species.set_id(0);
    auto ionisation_reaction = ElectronImpactIonisation<FixedRateData, FixedRateData>(
        A_particle_group->sycl_target, iz_rate_data, iz_rate_data, main_species,
        electron_species);

    // Recombination reaction
    // ------------------------------------------------------------------------------
    // Create Maxwellian distribution of recombination markers according to fluid
    // properties Add them to a new particle group representing the markers
    //

    // Options and constants

    // Make new particle group just for the markers
    auto marker_group =
        std::make_shared<ParticleGroup>(domain, particle_spec, sycl_target);

    // Give particle group initial kinetic values (positions and velocities)
    // Numerical settings: weight, stdev, species ID
    ParticleSet maxwellian_markers = uniform_cellwise_maxwellian<ndim>(
        sycl_target, neso_mesh, particle_spec, rec_markers_per_cell, 1.0, 0.5, -1);

    marker_group->add_particles_local(maxwellian_markers);

    // Give particle group initial fluid values: markers will contain background
    // plasma properties
    // From demo app "set_init_fluid_values"
    particle_loop(
        "set init fluid values", marker_group,
        [=](auto n, auto T, auto ne, auto Te, auto speed) {
          n.at(0) = background_ion_density;
          ne.at(0) = background_ion_density;
          T.at(0) = background_ion_temperature;
          Te.at(0) = background_ion_temperature;

          for (int i = 0; i < ndim; i++) {
            speed.at(i) = V_background[static_cast<std::size_t>(i)];
          }
        },
        Access::write(Sym<REAL>("FLUID_DENSITY")),
        Access::write(Sym<REAL>("FLUID_TEMPERATURE")),
        Access::write(Sym<REAL>("ELECTRON_DENSITY")),
        Access::write(Sym<REAL>("ELECTRON_TEMPERATURE")),
        Access::write(Sym<REAL>("FLUID_FLOW_SPEED")))
        ->execute();

    // Calculate marker weights

    // Add particle property: number of particles in the local cell
    // From demo app: "distribute_n_part_cell"
    for (int ic = 0; ic < marker_group->domain->mesh->get_cell_count(); ic++) {
      INT n_part_cell = marker_group->get_npart_cell(ic);
      particle_loop(
          "Update N_CELL prop", marker_group,
          [=](auto n_cell_prop) { n_cell_prop.at(0) = n_part_cell; },
          Access::write(Sym<INT>("N_CELL")))
          ->execute(ic);
    }

    const INT num_cells = marker_group->domain->mesh->get_cell_count();

    // Calculate weight for each marker particle
    // based on FLUID_DENSITY and N_CELL properties contained in same particle.
    // TODO: implement normalisation. dens_norm is currently 1.
    for (int ic = 0; ic < num_cells; ic++) {
      REAL V_cell = neso_mesh->dmh->get_cell_volume(ic);
      particle_loop(
          "Update weight of ions", marker_group,
          [=](auto n_cell_prop, auto ion_dens_prop, auto weight_prop) {
            const BoutReal n_cell = static_cast<BoutReal>(n_cell_prop.at(0));
            auto updated_weight = (ion_dens_prop.at(0) * Nnorm * V_cell) / (N_w * n_cell);
            weight_prop.at(0) = updated_weight;
          },
          Access::read(Sym<INT>("N_CELL")), Access::read(Sym<REAL>("FLUID_DENSITY")),
          Access::write(Sym<REAL>("WEIGHT")))
          ->execute(ic);
    }

    // Define marker species and reaction rates
    auto recomb_species = Species("ION", 1.0, 0.0, -1); // TODO: better as marker_species
    auto recomb_data = FixedRateData(rec_rate);
    auto recomb_energy_data = FixedRateData(rec_rate); // TODO: make this separate

    // This sampler will calculate marker momentum from fluid plasma conditions
    // TODO: Do I need a separate rng kernel?
    auto constant_rate_cross_section = ConstantRateCrossSection(1.0);

    auto recomb_data_calc_sampler =
        FilteredMaxwellianSampler<2, decltype(constant_rate_cross_section)>(
            1 / (recomb_species.get_mass() * Tnorm), constant_rate_cross_section,
            rng_kernel);

    // Container for objects allowing calculation of parameters within
    // the recombination kernel: sampled velocity and the radiation
    // energy loss source. Must be in this order.
    auto recomb_data_calc_obj =
        DataCalculator<decltype(recomb_energy_data), decltype(recomb_data_calc_sampler)>(
            recomb_energy_data, recomb_data_calc_sampler);

    BoutReal normalised_potential_energy = 1.0; // TODO: units
    auto recomb_reaction_kernel = RecombReactionKernels<2>(
        recomb_species, electron_species, normalised_potential_energy);

    // Set neutrals to be products of recombination
    const int out_state = static_cast<int>(main_species.get_id());
    std::array<int, 1> recomb_out_states = {out_state};

    // Create reaction object
    auto recomb_reaction =
        LinearReactionBase<1, decltype(recomb_data), decltype(recomb_reaction_kernel),
                           decltype(recomb_data_calc_obj)>(
            sycl_target, static_cast<int>(recomb_species.get_id()), recomb_out_states,
            recomb_data, recomb_reaction_kernel, recomb_data_calc_obj);

    // Wrappers & controllers
    // ------------------------------------------------------------------------------

    VantageSourceManager source_manager(neso_mesh, bout_mesh, units);

    const REAL remove_threshold = options["remove_threshold"].withDefault(1.0e-10);
    const REAL merge_threshold = options["merge_threshold"].withDefault(1.0e-2);

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

    // Ionisation reaction transforms and controller
    auto accumulator_transform_iz = std::make_shared<CellwiseAccumulator<REAL>>(
        A_particle_group, std::vector<std::string>{"ION_SOURCE_DENSITY"});

    auto accumulator_real_transform_wrapper = std::make_shared<TransformationWrapper>(
        std::dynamic_pointer_cast<TransformationStrategy>(accumulator_transform_iz));

    auto ion_source_density_zeroer =
        make_transformation_strategy<ParticleDatZeroer<REAL>>(
            std::vector<std::string>{"ION_SOURCE_DENSITY"});

    std::vector<std::shared_ptr<TransformationWrapper>> child_transforms =
        std::vector{merge_wrapper, remove_wrapper};
    std::vector<std::shared_ptr<TransformationWrapper>> parent_transforms_iz =
        std::vector{accumulator_real_transform_wrapper, merge_wrapper, remove_wrapper};

    auto reaction_controller = ReactionController(parent_transforms_iz, child_transforms);
    reaction_controller.add_reaction(
        std::make_shared<decltype(ionisation_reaction)>(ionisation_reaction));

    source_manager.add_source("Siz", "ION_SOURCE_DENSITY", accumulator_transform_iz,
                              A_particle_group, ion_source_density_zeroer);

    // Recombination transforms and controller

    auto accumulator_transform_rec = std::make_shared<CellwiseAccumulator<REAL>>(
        marker_group, std::vector<std::string>{"ION_SOURCE_DENSITY"});

    auto recomb_accumulator_transform_wrapper = std::make_shared<TransformationWrapper>(
        std::dynamic_pointer_cast<TransformationStrategy>(accumulator_transform_rec));

    std::vector<std::shared_ptr<TransformationWrapper>> parent_transforms_rec =
        std::vector{accumulator_real_transform_wrapper, merge_wrapper, remove_wrapper};

    auto recombination_controller =
        ReactionController(parent_transforms_rec, child_transforms);

    recombination_controller.add_reaction(
        std::make_shared<decltype(recomb_reaction)>(recomb_reaction));

    source_manager.add_source("Srec", "ION_SOURCE_DENSITY", accumulator_transform_rec,
                              marker_group, ion_source_density_zeroer);

    // Boundary handling
    // ------------------------------------------------------------------------------

    // Create the boundary interaction objects
    std::map<PetscInt, std::vector<PetscInt>> boundary_groups;
    // boundary_groups[1] = {100, 200};
    boundary_groups[1] = {100};

    b2d = std::make_shared<PetscInterface::BoundaryInteraction2D>(sycl_target, neso_mesh,
                                                                  boundary_groups);
    auto reflection = std::make_shared<BoundaryReflection>(ndim, 1.0e-10);

    auto lambda_apply_boundary_conditions = [&](auto aa) {
      auto sub_groups = b2d->post_integration(aa);
      for (auto& gx : sub_groups) {
        reflection->execute(gx.second, Sym<REAL>("POSITION"), Sym<REAL>("VELOCITY"),
                            Sym<REAL>("TSP"), b2d->previous_position_sym);
      }
    };

    // Advection & rest of code
    // ------------------------------------------------------------------------------

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
      const int size = static_cast<int>(get_npart_global(aa));
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
    h5part = std::make_shared<H5Part>(
        make_output_path("particle_trajectories.h5part", alloptions), A_particle_group,
        Sym<REAL>("POSITION"), Sym<REAL>("VELOCITY"));

    // allocate buffer vector for scalar projection/evaluation of NESO-Particles
    // properties
    std::vector<REAL> h_project1(static_cast<size_t>(num_cells_owned));
    // set weights from a Field2D from BOUT
    set_initial_particle_weights(initial_neutral_density, dg0, A_particle_group,
                                 neso_mesh, h_project1);

    // Calculate neutral density and sources for initial condition
    calculate_neutral_density_in_place(neutral_density, dg0, A_particle_group,
                                       h_project1);
    source_manager.update_all_sources(dt);

    // diagnose the initial condition
    std::string particle_data_filename =
        make_output_path(fmt::format("BOUT.dmp.vantage.{}.nc", mpi_rank), alloptions);
    Options bout_output_data = initialise_diagnostics(
        bout_mesh, neutral_density, ion_density, neso_mesh, particle_data_filename);
    // mass for conservation check
    Field2D total_density = neutral_density + ion_density;
    double total_mass_initial = calculate_total_mass(total_density, neso_mesh);
    // begin timestepping
    for (int stepx = 0; stepx < nsteps; stepx++) {
      // nprint("step:", stepx);
      output << "step:" << std::to_string(stepx) << std::endl;
      particle_time += dt;
      A_particle_group->hybrid_move();
      A_particle_group->cell_move();
      lambda_apply_timestep(static_particle_sub_group(A_particle_group));
      // apply reactions
      reaction_controller.apply(A_particle_group, dt, ControllerMode::standard_mode);
      recombination_controller.apply(marker_group, dt, A_particle_group);
      // uncomment to write a trajectory
      h5part->write();

      calculate_neutral_density_in_place(neutral_density, dg0, A_particle_group,
                                         h_project1);
      source_manager.update_all_sources(dt);
      Field2D Siz = source_manager.get_data("Siz");
      Field2D Srec = source_manager.get_data("Srec");

      // "Solve" density
      // Sources are in normalised m^-3 s^-1, so need to multiply by dt
      ion_density += (Siz + Srec) * dt;

      // diagnose timestep stepx
      update_diagnostics(neutral_density, ion_density, neso_mesh, bout_output_data,
                         particle_data_filename, particle_time);
    }
    // uncomment to write a trajectory

    // mass for conservation check
    total_density = neutral_density + ion_density;
    double total_mass_final = calculate_total_mass(total_density, neso_mesh);
    if (options["test_mass_conservation"].withDefault(true)) {
      check_mass_conservation(total_mass_final, total_mass_initial);
    }
  }
}

void Vantage::transform_impl(GuardedOptions& UNUSED(state)) {}

void Vantage::finally(const Options& UNUSED(state)) {}

void Vantage::outputVars(Options& UNUSED(state)) {}

// Destructor to handle VANTAGE related cleanup
Vantage::~Vantage() {
  h5part->close();          // Particle trajectory output
  b2d->free();              // NESO-Particles boundary interaction object
  neso_mesh->free();        // DMPlex interface
  PETSCCHK(DMDestroy(&dm)); // DMPlex mesh
  sycl_target->free();      // VANTAGE compute target
}

#endif
