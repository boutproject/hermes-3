
#pragma once
#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <functional>

#include "field3d.hxx"
#include "bout/coordinates.hxx"

/// Get the first argument from a parameter pack
template <typename Head, typename... Tail>
auto firstArg(const Head &head, Tail... ) {
  return head;
}

/// Take a function of BoutReals, and a region. Return
/// a function which takes fields (e.g. Field2D, Field3D),
/// and for every cell in the region evaluates the function
/// at quadrature points with weights.
/// These weights sum to 1, resulting in volume averaged values.
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
template <typename Function, typename RegionType>
auto cellAverage(Function func, RegionType region) {
  return [&](const auto &... args) {
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
      result[i] =
          4. / 6 * func((args[i])...) +
          (Ji + J[ym]) / (12. * Ji) * func((0.5 * (args[i] + args[ym]))...) +
          (Ji + J[yp]) / (12. * Ji) * func((0.5 * (args[i] + args[yp]))...);
    }
    return result;
  };
}

#endif // INTEGRATE_H
