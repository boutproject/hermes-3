#include "../include/component_scheduler.hxx"

#include "utils.hxx" // for trim, strsplit

#include "../include/ionisation.hxx"

ComponentScheduler::ComponentScheduler(Options &options,
                                       Solver *solver) {
  
  std::string component_names = options["components"]
                                    .doc("Components in order of execution")
                                    .as<std::string>();

  // For now split on ','. Something like "->" might be better
  for (const auto &name : strsplit(component_names, ',')) {
    // Ignore brackets, to allow these to be used to span lines.
    // In future brackets may be useful for complex scheduling
    auto name_trimmed = trim(name, " \t\r()");
    if (name_trimmed.empty()) {
      continue;
    }
    components.push_back(Component::create(name_trimmed,
                                           Options::root(),
                                           solver));
  }
}

std::unique_ptr<ComponentScheduler> ComponentScheduler::create(Options &options,
                                                               Solver *solver) {
  return std::make_unique<ComponentScheduler>(options, solver);
}


void ComponentScheduler::transform(Options &state) {
  // Run through each component
  for(auto &component : components) {
    component->transform(state);
  }
  // Enable components to update themselves based on the final state
  for(auto &component : components) {
    component->finally(state);
  }
}

void ComponentScheduler::annotate(Options &state) {
  // Run through each component
  for(auto &component : components) {
    component->annotate(state);
  }
}

void ComponentScheduler::precon(const Options &state, BoutReal gamma) {
  for(auto &component : components) {
    component->precon(state, gamma);
  }
}
