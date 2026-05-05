#pragma once
#include <string>
#ifndef VANTAGE_DMPLEX_H
#define VANTAGE_DMPLEX_H

#include "bout/bout.hxx"
#include "bout/bout_types.hxx"
#include "bout/field2d.hxx"
#include "bout/output.hxx"
#include "bout/petsclib.hxx"
#include <bout/field_factory.hxx>
#include <cmath>
#include <fmt/core.h>
#include <neso_particles.hpp>
#include <neso_particles/compute_target.hpp>
#include <neso_particles/containers/cell_data.hpp>
#include <neso_particles/external_interfaces/petsc/petsc_interface.hpp>
#include <neso_particles/typedefs.hpp>
#include <netcdf>
#include <petscsystypes.h>
#include <petscviewerhdf5.h>

#ifndef NESO_PARTICLES_PETSC
static_assert(false, "NESO-Particles was installed without PETSc support.");
#else

using namespace NESO::Particles;

void collect_unique_points(std::vector<double>& global_Z_vertices_buffer,
                           std::vector<double>& global_R_vertices_buffer, int& N_unique,
                           const double& tolerance,
                           std::vector<double>& global_Z_hypnotoad_vertices,
                           std::vector<double>& global_R_hypnotoad_vertices);

void RZ_to_ivertex_vector(Field2D& ivertex_corners,
                          std::vector<double>& global_Z_vertices,
                          std::vector<double>& global_R_vertices, const double& tolerance,
                          Mesh*& bout_mesh, Field2D& Rxy_corners, Field2D& Zxy_corners);

void load_vertex_information_from_netcdf(int& Nvertex,
                                         std::vector<double>& global_vertex_R,
                                         std::vector<double>& global_vertex_Z);

std::vector<PetscInt> cells_definition_from_RZ_ivertex(
    std::vector<PetscInt>& cells, Mesh*& bout_mesh, Field2D& Rxy_lower_left_corners,
    Field2D& Rxy_lower_right_corners, Field2D& Rxy_upper_right_corners,
    Field2D& Rxy_upper_left_corners, Field2D& Zxy_lower_left_corners,
    Field2D& Zxy_lower_right_corners, Field2D& Zxy_upper_right_corners,
    Field2D& Zxy_upper_left_corners, Field2D& ivertex_lower_left_corners,
    Field2D& ivertex_lower_right_corners, Field2D& ivertex_upper_right_corners,
    Field2D& ivertex_upper_left_corners);

DM create_dmplex_from_Bout_mesh(Mesh* bout_mesh, Options& mesh_options,
                                std::shared_ptr<SYCLTarget> sycl_target,
                                std::string dmplex_h5_filename);

#endif 

#endif // VANTAGE_DMPLEX_H
