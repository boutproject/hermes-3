
#pragma once

#ifndef COMPONENT_SCHEDULER_H
#define COMPONENT_SCHEDULER_H

#include "component.hxx"

#include <vector>
#include <memory>

/// Creates and schedules model components
///
/// Currently only one implementation, but in future alternative scheduler
/// types could be created. There is therefore a static create function
/// which in future could switch between types.
/// 
class ComponentScheduler {
public:
  ComponentScheduler(Options &scheduler_options, Options &component_options,
                     Solver *solver);

  /// Inputs
  ///   scheduler_options    Configuration of the scheduler
  ///     - components       Comma-separated list of component names
  ///
  ///   component_options    Configuration of the components.
  ///     - <name>
  ///       - ...
  ///   
  ///   solver               Used for time-dependent components
  ///                        to evolve quantities
  /// 
  static std::unique_ptr<ComponentScheduler> create(Options &scheduler_options,
                                                    Options &component_options,
                                                    Solver *solver);
  
  /// Run the scheduler, modifying the state
  void transform(Options &state);

  /// Add metadata, extra outputs
  void annotate(Options &state);

  /// Preconditioning
  void precon(const Options &state, BoutReal gamma);
private:
  /// The components to be executed in order
  std::vector<std::unique_ptr<Component>> components;
};

#endif // COMPONENT_SCHEDULER_H
