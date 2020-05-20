#pragma once
#ifndef ZERO_CURRENT
#define ZERO_CURRENT

#include "component.hxx"

/// Balance the parallel electron pressure gradient against
/// the electric field. Use this electric field to calculate
/// a force on the other species
///
///   E = -∇p_e / n_e
///
struct ZeroCurrent : public Component {
  ZeroCurrent(std::string, Options&, Solver*) {}
  
  /// Required inputs
  /// - species
  ///   - e
  ///     - pressure
  ///     - density
  ///     Asserts that charge = -1
  ///
  /// Sets in the input
  /// - species
  ///   - <all except e>   if both density and charge are set
  ///     - momentum_source
  /// 
  void transform(Options &state) override;
};

namespace {
RegisterComponent<ZeroCurrent> registercomponentzerocurrent("zero_current");
}

#endif // ZERO_CURRENT
