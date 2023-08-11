
#pragma once
#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <functional>

#include <bout/field3d.hxx>
#include <bout/coordinates.hxx>
#include <bout/fv_ops.hxx>

#include "../include/hermes_build_config.hxx"

/// Get the first argument from a parameter pack
template <typename Head, typename... Tail>
auto firstArg(const Head &head, Tail... ) {
  return head;
}

/// Return the value at the left of a cell,
/// given cell centre values at this cell and two neighbours
template <typename CellEdges>
BoutReal cellLeft(BoutReal c, BoutReal m, BoutReal p) {
  CellEdges cellboundary;
  FV::Stencil1D s {.c = c, .m = m, .p = p};
  cellboundary(s);
  return s.L;
}

/// Return the value at the right of a cell,
/// given cell centre values at this cell and two neighbours
template <typename CellEdges>
BoutReal cellRight(BoutReal c, BoutReal m, BoutReal p) {
  CellEdges cellboundary;
  FV::Stencil1D s {.c = c, .m = m, .p = p};
  cellboundary(s);
  return s.R;
}

/// Take a function of BoutReals, and a region. Return
/// a function which takes fields (e.g. Field2D, Field3D),
/// and for every cell in the region evaluates the function
/// at quadrature points with weights.
/// These weights sum to 1, resulting in volume averaged values.
///
/// Uses a limiter to calculate values at cell edges. This is
/// needed so that as Ne goes to zero in a cell then atomic
/// rates also go to zero.
///
/// Example
///   Field3D Ne = ..., Te = ...;
/// 
///   Field3D result = cellAverage(
///          [](BoutReal Ne, BoutReal Te) {return Ne*Te;} // The function to evaluate
///          Ne.getRegion("RGN_NOBNDRY")  // The region to iterate over
///          )(Ne, Te);                   // The input fields
///
/// Note that the order of the arguments to the lambda function
/// is the same as the input fields.
///
template <typename CellEdges = hermes::Limiter, typename Function, typename RegionType>
auto cellAverage(Function func, const RegionType &region) {
  // Note: Capture by value or func and region go out of scope
  return [=](const auto &... args) { 
    // Use the first argument to set the result mesh etc.
    Field3D result{emptyFrom(firstArg(args...))};
    result.allocate();

    // Get the coordinate Jacobian
    auto J = result.getCoordinates()->J;
    BOUT_FOR(i, region) {
      // Offset indices
      auto yp = i.yp();
      auto ym = i.ym();
      auto Ji = J[i];

      // Integrate in Y using Simpson's rule
      // Using limiter to calculate cell edge values
      result[i] =
        4. / 6 * func((args[i])...) +
        (Ji + J[ym]) / (12. * Ji) * func(cellLeft<CellEdges>(args[i], args[ym], args[yp])...) +
        (Ji + J[yp]) / (12. * Ji) * func(cellRight<CellEdges>(args[i], args[ym], args[yp])...);
    }
    return result;
  };
}

#endif // INTEGRATE_H
