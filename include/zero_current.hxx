#pragma once
#ifndef ZERO_CURRENT
#define ZERO_CURRENT

#include "component.hxx"

/// Balance the parallel electron pressure gradient against
/// the electric field. Use this electric field to calculate
/// a force on the other species
///
///   E = (-∇p_e + F) / n_e
///
/// where F is the momentum source for the electrons
///
struct ZeroCurrent : public Component {
  ZeroCurrent(std::string, Options&, Solver*) {}
  
  /// Required inputs
  /// - species
  ///   - e
  ///     - pressure
  ///     - density
  ///     - momentum_source [optional]
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
