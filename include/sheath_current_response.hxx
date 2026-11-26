#pragma once

#include <bout/field2d.hxx>

/// Toroidal mean of the electron sheath-face momentum derivative with respect
/// to the face potential, stored in each adjacent interior cell.
struct SheathCurrentResponse {
  const Field2D& enableCurrentResponse() {
    current_response = 0.0;
    return current_response;
  }

protected:
  Field2D current_response; ///< Allocated only when requested by the n=0 solve
};
