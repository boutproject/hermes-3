#pragma once
#include "bout/bout.hxx"
#include <neso_particles.hpp>
#include <reactions/reactions.hpp>

using namespace NESO::Particles;
using namespace VANTAGE::Reactions;

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