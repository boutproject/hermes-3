
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
  ComponentScheduler(Options &options, const MeshMap &meshes);

  static std::unique_ptr<ComponentScheduler> create(Options &options,
                                                    const MeshMap &meshes);
  
  /// Run the scheduler, modifying the state
  void transform(Options &state);
  
private:
  /// The components to be executed in order
  std::vector<std::unique_ptr<Component>> components;
};

#endif // COMPONENT_SCHEDULER_H
