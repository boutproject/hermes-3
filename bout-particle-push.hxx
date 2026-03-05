#pragma once
#include "bout/bout.hxx"
#include <neso_particles.hpp>
#include <neso_rng_toolkit.hpp>
#include <reactions/reactions.hpp>

using namespace NESO::Particles;
using namespace VANTAGE::Reactions;

/**
 * @brief Function to calculate cell volumes.
 * 
 * @param sycl_target SYCLTargetSharedPtr to use for communication.
 * @param mesh object.
 * 
 */


inline auto calc_V_tot_local(SYCLTargetSharedPtr sycl_target,
                             std::shared_ptr<PetscInterface::DMPlexInterface> mesh,
                             std::map<std::string, double> &norms) {
  auto num_cells = mesh->get_cell_count();

  std::vector<REAL> V_cells;
  REAL V_tot_local = 0.0;
  REAL V_tot_global;
  for (int cellx = 0; cellx < num_cells; cellx++) {
    // No normalisation needed unlike in the demo app, our grid is in physical units
    const REAL V_cell = mesh->dmh->get_cell_volume(cellx) * norms["length"] * norms["length"];
    V_cells.push_back(V_cell);
    V_tot_local += V_cell;
  }

  MPI_Allreduce(&V_tot_local, &V_tot_global, 1, MPI_DOUBLE, MPI_SUM,
                sycl_target->comm_pair.comm_parent);

  const int rank = sycl_target->comm_pair.rank_parent;

  if (rank == 0) {
    std::cout << "Length_norm: " << norms["length"] << std::endl;
    std::cout << "V_tot_global: " << V_tot_global << std::endl;
  }

  return std::tuple(V_tot_local, V_tot_global, V_cells);
}


/**
 * @brief Function to calculate particle positions and velocities from a Maxwellian.
 * 
 * @param sycl_target SYCLTargetSharedPtr to use for communication.
 * @param mesh object.
 * @param particle_spec ParticleSpec to use for the returned ParticleSet.
 * @param npart_per_cell Number of particles to create per cell.
 * @param weight Weight to assign to each particle.
 * @param std_dev Standard deviation of the Maxwellian distribution to sample
 * velocities from.
 * @param species_id Integer to assign to the "INTERNAL_STATE" property of each
 * particle.
 */

template <size_t ndim>
inline ParticleSet uniform_cellwise_maxwellian(
    SYCLTargetSharedPtr sycl_target, std::shared_ptr<PetscInterface::DMPlexInterface> mesh,
    const ParticleSpec &particle_spec, const INT &npart_per_cell,
    const REAL &weight, const REAL &std_dev, const INT &species_id) {

  const int rank = sycl_target->comm_pair.rank_parent;
  const int size = sycl_target->comm_pair.size_parent;

  std::mt19937 rng_pos(52234234 + rank);
  std::mt19937 rng_vel(52234231 + rank);

  std::vector<std::vector<double>> positions;
  std::vector<int> cell_ids;
  PetscInterface::uniform_within_dmplex_cells(mesh, npart_per_cell, positions, cell_ids,
                                 &rng_pos);

  const int N = cell_ids.size();

  auto velocities =
      NESO::Particles::normal_distribution(N, ndim, 0.0, std_dev, rng_vel);

  ParticleSet maxwellian(N, particle_spec);

  for (int px = 0; px < N; px++) {
    for (int dimx = 0; dimx < ndim; dimx++) {
      maxwellian[Sym<REAL>("POSITION")][px][dimx] = positions.at(dimx).at(px);
      maxwellian[Sym<REAL>("VELOCITY")][px][dimx] = velocities.at(dimx).at(px);
    }
    maxwellian[Sym<INT>("CELL_ID")][px][0] = cell_ids.at(px);
    maxwellian[Sym<INT>("ID")][px][0] = px;
    maxwellian[Sym<REAL>("WEIGHT")][px][0] = weight;
    maxwellian[Sym<INT>("INTERNAL_STATE")][px][0] = species_id;
  }

  return maxwellian;
}

/**
 * @brief create an RNG kernel to use for sampling velocity distributions
 * in recombination and charge exchange. 
 */

inline auto get_uniform_rng_kernel(SYCLTargetSharedPtr sycl_target,
                                   std::size_t n_samples,
                                   std::uint64_t root_seed = 141351) {

  const int rank = sycl_target->comm_pair.rank_parent;

  std::uint64_t seed = NESO::RNGToolkit::create_seeds(
      sycl_target->comm_pair.size_parent, rank, root_seed);

  auto rng_normal = NESO::RNGToolkit::create_rng<REAL>(
      NESO::RNGToolkit::Distribution::Uniform<REAL>{
          NESO::RNGToolkit::Distribution::next_value(0.0), 1.0},
      seed, sycl_target->device, sycl_target->device_index);

  // Create an interface between NESO-RNG-Toolkit and NESO-Particles KernelRNG
  auto rng_interface =
      make_rng_generation_function<GenericDeviceRNGGenerationFunction, REAL>(
          [=](REAL *d_ptr, const std::size_t num_samples) -> int {
            return rng_normal->get_samples(d_ptr, num_samples);
          });

  auto rng_kernel =
      host_atomic_block_kernel_rng<REAL>(rng_interface, n_samples);

  return rng_kernel;
}