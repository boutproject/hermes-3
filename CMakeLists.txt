cmake_minimum_required(VERSION 3.9...3.12)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
  cmake_policy(VERSION 3.12)
endif()

project(hermes-3 LANGUAGES CXX C)

# Extra CMake scripts in cmake/ subdirectory
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

option(HERMES_UPDATE_GIT_SUBMODULE "Check submodules are up-to-date during build" ON)
# Adapted from https://cliutils.gitlab.io/modern-cmake/chapters/projects/submodule.html
# Update submodules as needed
function(hermes_update_submodules)
  if(NOT HERMES_UPDATE_GIT_SUBMODULE)
    return()
  endif()
  find_package(Git QUIET)
  if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    message(STATUS "Submodule update")
    execute_process(COMMAND ${GIT_EXECUTABLE} -c submodule.recurse=false submodule update --init --recursive
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()
  endif()
endfunction()
hermes_update_submodules()

# Get the Git revision
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC HERMES_REVISION)
if(HERMES_REVISION STREQUAL "GITDIR-NOTFOUND")
  set(HERMES_REVISION "Unknown")
endif()
message(STATUS "Git revision: ${HERMES_REVISION}")

# BOUT++ is a dependency
option(HERMES_BUILD_BOUT "Build BOUT++ in external/BOUT-dev" ON)
if(HERMES_BUILD_BOUT)
  set(BOUT_BUILD_EXAMPLES, OFF) # Don't create example makefiles
  set(HERMES_BOUT_SRC "external/BOUT-dev" CACHE STRING "Directory of BOUT++ dir")
  add_subdirectory(${HERMES_BOUT_SRC} ${CMAKE_CURRENT_BINARY_DIR}/external/BOUT-dev)
else()
  find_package(bout++ REQUIRED)
endif()

set(HERMES_SOURCES
    src/classical_diffusion.cxx
    src/component.cxx
    src/component_scheduler.cxx
    src/ionisation.cxx
    src/radiation.cxx
    src/div_ops.cxx
    src/neutral_mixed.cxx
    src/full_velocity.cxx
    src/electromagnetic.cxx
    src/electron_force_balance.cxx
    src/electron_viscosity.cxx
    src/evolve_density.cxx
    src/evolve_energy.cxx
    src/evolve_pressure.cxx
    src/evolve_momentum.cxx
    src/isothermal.cxx
    src/quasineutral.cxx
    src/diamagnetic_drift.cxx
    src/recalculate_metric.cxx
    src/relax_potential.cxx
    src/sheath_closure.cxx
    src/sheath_boundary.cxx
    src/sheath_boundary_simple.cxx
    src/sheath_boundary_insulating.cxx
    src/snb_conduction.cxx
    src/fixed_fraction_ions.cxx
    src/sound_speed.cxx
    src/zero_current.cxx
    src/collisions.cxx
    src/anomalous_diffusion.cxx
    src/binormal_stpm.cxx
    src/recycling.cxx
    src/amjuel_hyd_ionisation.cxx
    src/amjuel_hyd_recombination.cxx
    src/amjuel_helium.cxx
    src/adas_reaction.cxx
    src/noflow_boundary.cxx
    src/neutral_parallel_diffusion.cxx
    src/neutral_boundary.cxx
    src/polarisation_drift.cxx
    src/solkit_neutral_parallel_diffusion.cxx
    src/hydrogen_charge_exchange.cxx
    src/solkit_hydrogen_charge_exchange.cxx
    src/upstream_density_feedback.cxx
    src/temperature_feedback.cxx
    src/detachment_controller.cxx
    src/thermal_force.cxx
    src/ion_viscosity.cxx
    src/transform.cxx
    src/vorticity.cxx
    include/adas_reaction.hxx
    include/adas_carbon.hxx
    include/adas_neon.hxx
    include/amjuel_helium.hxx
    include/amjuel_hyd_ionisation.hxx
    include/amjuel_hyd_recombination.hxx
    include/amjuel_reaction.hxx
    include/anomalous_diffusion.hxx
    include/classical_diffusion.hxx
    include/binormal_stpm.hxx
    include/collisions.hxx
    include/component.hxx
    include/component_scheduler.hxx
    include/diamagnetic_drift.hxx
    include/div_ops.hxx
    include/electromagnetic.hxx
    include/electron_force_balance.hxx
    include/electron_viscosity.hxx
    include/evolve_density.hxx
    include/evolve_energy.hxx
    include/evolve_momentum.hxx
    include/evolve_pressure.hxx
    include/fixed_density.hxx
    include/fixed_fraction_ions.hxx
    include/fixed_velocity.hxx
    include/full_velocity.hxx
    include/hermes_utils.hxx
    include/hydrogen_charge_exchange.hxx
    include/solkit_hydrogen_charge_exchange.hxx
    include/integrate.hxx
    include/ion_viscosity.hxx
    include/ionisation.hxx
    include/neutral_boundary.hxx
    include/neutral_mixed.hxx
    include/neutral_parallel_diffusion.hxx
    include/solkit_neutral_parallel_diffusion.hxx
    include/noflow_boundary.hxx
    include/polarisation_drift.hxx
    include/quasineutral.hxx
    include/radiation.hxx
    include/recalculate_metric.hxx
    include/recycling.hxx
    include/relax_potential.hxx
    include/sheath_boundary.hxx
    include/sheath_boundary_simple.hxx
    include/sheath_boundary_insulating.hxx
    include/simple_conduction.hxx
    include/snb_conduction.hxx
    include/sheath_closure.hxx
    include/sound_speed.hxx
    include/thermal_force.hxx
    include/upstream_density_feedback.hxx
    include/temperature_feedback.hxx
    include/detachment_controller.hxx
    include/vorticity.hxx
    include/zero_current.hxx
    include/transform.hxx
    include/fixed_fraction_radiation.hxx
    include/simple_pump.hxx
    )

# Compile Hermes sources into a library, which can be used for unit tests
add_library(hermes-3-lib ${HERMES_SOURCES})

target_link_libraries(hermes-3-lib PRIVATE bout++::bout++)

target_include_directories(hermes-3-lib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>
  )

# The main executable target
add_executable(hermes-3
               hermes-3.cxx
               ${CMAKE_CURRENT_BINARY_DIR}/include/revision.hxx)

target_link_libraries(hermes-3 PRIVATE bout++::bout++ hermes-3-lib)

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

# Once built, copy the data and test directories
option(HERMES_COPY_EXAMPLES_TO_BUILD "Copy the examples to the build directory" ON)
if(HERMES_COPY_EXAMPLES_TO_BUILD)
  add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E copy_directory
                    ${CMAKE_SOURCE_DIR}/examples $<TARGET_FILE_DIR:${PROJECT_NAME}>/examples)
endif()

option(HERMES_COPY_JSON_DATABASE_TO_BUILD "Copy the json database to the build directory" ON)
if(HERMES_COPY_JSON_DATABASE_TO_BUILD)
  add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E copy_directory
                    ${CMAKE_SOURCE_DIR}/json_database $<TARGET_FILE_DIR:${PROJECT_NAME}>/json_database)
endif()

option(HERMES_COPY_TESTS_TO_BUILD "Copy the tests to the build directory" ON)
if(HERMES_COPY_TESTS_TO_BUILD)
  add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E copy_directory
                    ${CMAKE_SOURCE_DIR}/tests $<TARGET_FILE_DIR:${PROJECT_NAME}>/tests)
endif()

# Tests
option(HERMES_TESTS "Build the tests" ON)
if(HERMES_TESTS)
  enable_testing()

  # Integrated tests
  function(hermes_add_integrated_test TESTNAME)
    add_test(NAME ${TESTNAME}
           WORKING_DIRECTORY tests/integrated/${TESTNAME}
           COMMAND runtest)
  endfunction()

  hermes_add_integrated_test(1D-fluid)
  hermes_add_integrated_test(1D-recycling)
  hermes_add_integrated_test(diffusion)
  hermes_add_integrated_test(evolve_density)
  hermes_add_integrated_test(neutral_mixed)
  hermes_add_integrated_test(vorticity)
  hermes_add_integrated_test(sod-shock)
  hermes_add_integrated_test(sod-shock-energy)
  hermes_add_integrated_test(drift-wave)
  hermes_add_integrated_test(alfven-wave)

  # Unit tests
  option(HERMES_UNIT_TESTS "Build the unit tests" ON)
  if(HERMES_UNIT_TESTS)
    if(NOT PACKAGE_TESTS)
      # Not running BOUT++ tests
      # Since BOUT++ includes googletest, adding it again causes errors
      set(CMAKE_CXX_STANDARD 14) # Needed to compile GoogleTest on Darwin
      add_subdirectory(external/googletest)
    endif()
    file(GLOB_RECURSE TEST_SOURCES tests/unit/test_*.cxx)

    add_executable(hermes_unit_tests
                   tests/unit/main.cxx
                   ${TEST_SOURCES})

    add_test(NAME hermes_unit_tests COMMAND hermes_unit_tests)
    target_link_libraries(hermes_unit_tests PRIVATE gtest bout++::bout++ hermes-3-lib)
  endif()
endif()

# Compile-time options

set(SLOPE_LIMITERS MinMod MC Upwind Superbee)
set(HERMES_SLOPE_LIMITER MC CACHE STRING "Set advection slope limiter")
set_property(CACHE HERMES_SLOPE_LIMITER PROPERTY STRINGS ${SLOPE_LIMITERS})
message(STATUS "Slope limiter: ${HERMES_SLOPE_LIMITER}")

# Generate the build config header

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include/hermes_build_config.hxx")
  # If we do in source builds, this is fine
  if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
    message(FATAL_ERROR "Generated hermes_build_config.hxx header already exists; please remove '${CMAKE_CURRENT_SOURCE_DIR}/include/hermes_build_config.hxx' before continuing")
  endif()
endif()

configure_file(include/hermes_build_config.hxx.in include/hermes_build_config.hxx)

# Installation
if (BUILD_SHARED_LIBS)
  install(TARGETS hermes-3 hermes-3-lib)
else()
  install(TARGETS hermes-3)
endif()
