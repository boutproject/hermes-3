#include <bout/boutcomm.hxx>
#include "gtest/gtest.h"

int main(int argc, char** argv) {
  // Newer BOUT++ initialises MPI lazily via `BoutComm::getComm()`, and will call
  // `MPI_Init(pargc, pargv)` when first needed. Make sure those pointers are set
  // so we don't pass nulls into MPI.
  BoutComm::setArgs(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
