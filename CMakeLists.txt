cmake_minimum_required(VERSION 3.9...3.12)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
  cmake_policy(VERSION 3.12)
endif()

project(hermes-3 LANGUAGES CXX)

# Extra CMake scripts in cmake/ subdirectory
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Update submodules
# Adapted from https://cliutils.gitlab.io/modern-cmake/chapters/projects/submodule.html
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  message(STATUS "Submodule update")
  execute_process(COMMAND ${GIT_EXECUTABLE} -c submodule.recurse=false submodule update --init --recursive
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    RESULT_VARIABLE GIT_SUBMOD_RESULT)
  if(NOT GIT_SUBMOD_RESULT EQUAL "0")
    message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
  endif()
endif()

# Get the Git revision
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC HERMES_REVISION)
if(HERMES_REVISION STREQUAL "GITDIR-NOTFOUND")
  set(HERMES_REVISION "Unknown")
endif()
message(STATUS "Git revision: ${HERMES_REVISION}")

# BOUT++ is a dependency
add_subdirectory(external/BOUT-dev)

add_executable(hermes-3
               hermes-3.cxx
               src/component.cxx
               src/component_scheduler.cxx
               src/ionisation.cxx
               src/radiation.cxx
               src/div_ops.cxx
               src/loadmetric.cxx
               src/neutral_mixed.cxx
               src/full_velocity.cxx
               src/evolve_density.cxx
               src/vorticity.cxx
               src/evolve_pressure.cxx
               src/evolve_momentum.cxx
               src/isothermal.cxx
               src/quasineutral.cxx
               src/diamagnetic_drift.cxx
               src/sheath_closure.cxx
               src/sheath_boundary.cxx
               src/fixed_fraction_ions.cxx
               src/sound_speed.cxx
               src/zero_current.cxx
               src/collisions.cxx
               src/anomalous_diffusion.cxx
               src/recycling.cxx
               src/amjuel_hyd_ionisation.cxx
               src/amjuel_helium.cxx
               src/adas_reaction.cxx
               src/noflow_boundary.cxx
               src/neutral_parallel_diffusion.cxx
               src/hydrogen_charge_exchange.cxx
               src/upstream_density_feedback.cxx
               src/thermal_force.cxx
               include/adas_reaction.hxx
               include/adas_neon.hxx
               include/amjuel_helium.hxx
               include/amjuel_hyd_ionisation.hxx
               include/amjuel_reaction.hxx
               include/anomalous_diffusion.hxx
               include/collisions.hxx
               include/component.hxx
               include/component_scheduler.hxx
               include/diamagnetic_drift.hxx
               include/div_ops.hxx
               include/evolve_density.hxx
               include/evolve_momentum.hxx
               include/evolve_pressure.hxx
               include/fixed_fraction_ions.hxx
               include/full_velocity.hxx
               include/hydrogen_charge_exchange.hxx
               include/integrate.hxx
               include/ionisation.hxx
               include/loadmetric.hxx
               include/neutral_mixed.hxx
               include/neutral_parallel_diffusion.hxx
               include/noflow_boundary.hxx
               include/quasineutral.hxx
               include/radiation.hxx
               include/recycling.hxx
               include/sheath_boundary.hxx
               include/sheath_closure.hxx
               include/sound_speed.hxx
               include/thermal_force.hxx
               include/upstream_density_feedback.hxx
               include/vorticity.hxx
               include/zero_current.hxx
               ${CMAKE_CURRENT_BINARY_DIR}/include/revision.hxx
               )

target_link_libraries(hermes-3 PRIVATE bout++::bout++)
target_include_directories(hermes-3 PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>
  )

# Build the file containing just the commit hash
# This will be rebuilt on every commit!
configure_file(
  "${PROJECT_SOURCE_DIR}/include/revision.hxx.in"
  "${PROJECT_BINARY_DIR}/include/revision.hxx")
