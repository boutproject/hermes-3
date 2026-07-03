#pragma once
#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <bout/bout_types.hxx>
#include <bout/coordinates.hxx>
#include <bout/field3d.hxx>
#include <bout/fv_ops.hxx>

#include "../include/hermes_build_config.hxx"

/// Get the first argument from a parameter pack
template <typename Head, typename... Tail>
auto firstArg(const Head& head, Tail...) {
  return head;
}

inline auto threePointStencil(BoutReal c, BoutReal m, BoutReal p) {
  // TODO(peter): We can remove this #ifdef guard after switching to C++20
#if __cpp_designated_initializers >= 201707L
  return FV::Stencil1D{
      .c = c, .m = m, .p = p, .mm = BoutNaN, .pp = BoutNaN, .L = BoutNaN, .R = BoutNaN};
#else
  return FV::Stencil1D{c, m, p, BoutNaN, BoutNaN, BoutNaN, BoutNaN};
#endif
}

/// Return the value at the left of a cell,
/// given cell centre values at this cell and two neighbours
template <typename CellEdges>
BoutReal cellLeft(BoutReal c, BoutReal m, BoutReal p) {
  CellEdges cellboundary;
  FV::Stencil1D s = threePointStencil(c, m, p);
  cellboundary(s);
  return s.L;
}

/// Return the value at the right of a cell,
/// given cell centre values at this cell and two neighbours
template <typename CellEdges>
BoutReal cellRight(BoutReal c, BoutReal m, BoutReal p) {
  CellEdges cellboundary;
  FV::Stencil1D s = threePointStencil(c, m, p);
  cellboundary(s);
  return s.R;
}

/// Take a function of BoutReals, a region and input fields
/// (e.g. Field2D, Field3D). For every cell in the region evaluate
/// the function at quadrature points with weights.
/// These weights sum to 1, resulting in volume averaged values,
/// written into the caller-supplied `result` field.
///
/// Reusing `result` between calls avoids constructing a new Field3D
/// on every call: on the first call it is initialised from the first
/// input field; subsequent calls reuse the existing data.
///
/// Uses a limiter to calculate values at cell edges. This is
/// needed so that as Ne goes to zero in a cell then atomic
/// rates also go to zero.
///
/// Note that the order of the arguments to the lambda function
/// is the same as the input fields.
template <typename CellEdges = hermes::Limiter, typename Function, typename RegionType,
          typename... Fields>
void cellAverageInto(Field3D& result, Function func, const RegionType& region,
                     const Fields&... args) {
  if (!result.isAllocated()) {
    // Use the first argument to set the result mesh etc.
    result = emptyFrom(firstArg(args...));
  }
  // Ensure the backing data is not shared with another field (copy-on-write)
  result.allocate();

  // Get the coordinate Jacobian
  auto J = result.getCoordinates()->J;
  BOUT_FOR(i, region) {
    // Offset indices
    auto yp = i.yp();
    auto ym = i.ym();
    auto Ji = J[i];
    const BoutReal inv_12_Ji = 1.0 / (12.0 * Ji);

    // Integrate in Y using Simpson's rule
    // Using limiter to calculate cell edge values
    result[i] = (2. / 3) * func((args[i])...)
                + (Ji + J[ym]) * inv_12_Ji
                      * func(cellLeft<CellEdges>(args[i], args[ym], args[yp])...)
                + (Ji + J[yp]) * inv_12_Ji
                      * func(cellRight<CellEdges>(args[i], args[ym], args[yp])...);
  }
}

/// Like cellAverageInto, but takes only the function and region, and
/// returns a function which takes the input fields and returns the
/// volume-averaged result as a new Field3D.
///
/// Example
///   Field3D Ne = ..., Te = ...;
///
///   Field3D result = cellAverage(
///          [](BoutReal Ne, BoutReal Te) {return Ne*Te;} // The function to evaluate
///          Ne.getRegion("RGN_NOBNDRY")  // The region to iterate over
///          )(Ne, Te);                   // The input fields
///
template <typename CellEdges = hermes::Limiter, typename Function, typename RegionType>
auto cellAverage(Function func, const RegionType& region) {
  // Note: Capture by value or func and region go out of scope
  return [=](const auto&... args) {
    Field3D result;
    cellAverageInto<CellEdges>(result, func, region, args...);
    return result;
  };
}

#endif // INTEGRATE_H
