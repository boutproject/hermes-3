
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
  ///  @param scheduler_options  Configuration of the scheduler
  ///                            Should contain "components", a comma-separated
  ///                            list of component names
  ///
  ///  @param component_options  Configuration of the components.
  ///     - <name>
  ///       - type = Component classes, ...
  ///                        If not provided, the type is the same as the name
  ///                        Multiple classes can be given, separated by commas.
  ///                        All classes will use the same Options section.
  ///       - ...  Options to control the component(s)
  ///   
  ///  @param solver         Used for time-dependent components
  ///                        to evolve quantities
  /// 
  static std::unique_ptr<ComponentScheduler> create(Options &scheduler_options,
                                                    Options &component_options,
                                                    Solver *solver);
  
  /// Run the scheduler, modifying the state.
  /// This calls all components' transform() methods, then
  /// all component's finally() methods.
  void transform(Options &state);

  /// Add metadata, extra outputs. This would typically
  /// be called only for writing to disk, rather than every internal
  /// timestep.
  void outputVars(Options &state);

  /// Add variables to restart files
  void restartVars(Options &state);

  /// Preconditioning
  void precon(const Options &state, BoutReal gamma);
private:
  /// The components to be executed in order
  std::vector<std::unique_ptr<Component>> components;
};

#endif // COMPONENT_SCHEDULER_H
