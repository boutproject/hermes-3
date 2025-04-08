set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")

option(HERMES_USE_SCOREP "Enable support for Score-P based instrumentation" OFF)
if (HERMES_USE_SCOREP)
  message(STATUS "Score-P support enabled. Please make sure you are calling CMake like so:

  SCOREP_WRAPPER=off cmake -DCMAKE_C_COMPILER=scorep-mpicc -DCMAKE_CXX_COMPILER=scorep-mpicxx <other CMake options>
")
endif()
set(HERMES_HAS_SCOREP ${HERMES_USE_SCOREP})