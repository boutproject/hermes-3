#include <memory>
#include <string>
#include <vector>

#include <bout/bout_types.hxx>
#include <bout/options.hxx>
#include <bout/utils.hxx> // for trim, strsplit

#include "../include/component.hxx"
#include "../include/component_scheduler.hxx"

/// Perform a depth-first topological sort, starting from `item`.
void topological_sort(const std::vector<std::set<size_t>>& dependencies, size_t item,
                      std::vector<size_t>& sorted, std::vector<bool>& processing,
                      std::vector<bool>& processed) {
  if (processed[item]) {
    return;
  }
  if (processing[item]) {
    throw BoutException("Circular dependency among components.");
  }
  processing[item] = true;

  for (const auto dep : dependencies[item]) {
    topological_sort(dependencies, dep, sorted, processing, processed);
  }
  processed[item] = true;
  sorted.push_back(item);
}

/// Consumes a list of components and returns a new one that has been
/// topolgically sorted to ensure variables are written and read in
/// the right order.
std::vector<std::unique_ptr<Component>>
sortComponents(std::vector<std::unique_ptr<Component>>&& components) {
  using Var = std::pair<std::string, Regions>;
  std::map<Var, std::set<size_t>> nonfinal_writes, final_writes;
  // Build up information on which components write each variable
  for (size_t i = 0; i < components.size(); i++) {
    const Permissions& permissions = components[i]->getPermissions();
    for (const auto& [name, regions] :
         permissions.getVariablesWithPermission(PermissionTypes::Write)) {
      for (const auto& [region, _] : Permissions::fundamental_regions) {
        if ((regions & region) == region) {
          nonfinal_writes[{name, region}].insert(i);
        }
      }
    }
    for (const auto& [name, regions] :
         permissions.getVariablesWithPermission(PermissionTypes::Final)) {
      for (const auto& [region, _] : Permissions::fundamental_regions) {
        if ((regions & region) == region) {
          final_writes[{name, region}].insert(i);
        }
      }
    }
  }

  std::vector<std::set<size_t>> component_dependencies(components.size());

  // Components which do a final write on a variable depend on all
  // components which do non-final writes on that variable
  for (const auto& [var, comp_indices] : final_writes) {
    for (size_t i : comp_indices) {
      const auto item = nonfinal_writes.find(var);
      if (item != nonfinal_writes.end()) {
        component_dependencies[i].merge(item->second);
      }
    }
  }

  // Work out which component(s) last write a variable before it may be read
  std::map<Var, std::set<size_t>> variable_writers = std::move(final_writes);
  variable_writers.merge(std::move(nonfinal_writes));

  for (size_t i = 0; i < components.size(); i++) {
    const Permissions& permissions = components[i]->getPermissions();
    // Create dependencies between components that read variables and those that write
    // them
    for (const auto& [name, regions] :
         permissions.getVariablesWithPermission(PermissionTypes::Read)) {
      for (const auto& [region, _] : Permissions::fundamental_regions) {
        if ((regions & region) == region) {
          const auto item = variable_writers.find({name, region});
          if (item == variable_writers.end()) {
            throw BoutException(
                "Required variable {} (in region {}) is not written by any component.",
                name, Permissions::fundamental_regions.at(region));
          }
          component_dependencies[i].merge(item->second);
        }
      }
    }

    // Create dependencies for ReadIfSet variables only if there exist another component
    // which sets them
    for (const auto& [name, regions] :
         permissions.getVariablesWithPermission(PermissionTypes::ReadIfSet)) {
      for (const auto& [region, _] : Permissions::fundamental_regions) {
        if ((regions & region) == region) {
          const auto item = variable_writers.find({name, region});
          if (item != variable_writers.end()) {
            component_dependencies[i].merge(item->second);
          }
        }
      }
    }
  }

  std::vector<bool> processing(components.size(), false),
      processed(components.size(), false);
  std::vector<std::size_t> order;

  for (size_t i = 0; i < components.size(); i++) {
    if (!processed[i]) {
      topological_sort(component_dependencies, i, order, processing, processed);
    }
  }

  // Create the result with components in the desired order
  std::vector<std::unique_ptr<Component>> result(components.size());
  for (size_t i = 0; i < components.size(); i++) {
    std::swap(result[i], components[order[i]]);
  }

  return result;
}

ComponentScheduler::ComponentScheduler(Options &scheduler_options,
                                       Options &component_options,
                                       Solver *solver) {

  const std::string component_names = scheduler_options["components"]
                                          .doc("Components in order of execution")
                                          .as<std::string>();

  std::vector<std::string> electrons;
  std::vector<std::string> neutrals;
  std::vector<std::string> positive_ions;
  std::vector<std::string> negative_ions;

  // For now split on ','. Something like "->" might be better
  for (const auto &name : strsplit(component_names, ',')) {
    // Ignore brackets, to allow these to be used to span lines.
    // In future brackets may be useful for complex scheduling

    auto name_trimmed = trim(name, " \t\r()");
    if (name_trimmed.empty()) {
      continue;
    }

    if (name_trimmed == "e" or name == "ebeam") {
      electrons.push_back(name_trimmed);
    }
    // FIXME: Would there be any spcies without AA? Is there any other
    // reliable way to identify what is a species?
    else if (component_options[name_trimmed].isSet("AA")) {
      if (component_options[name_trimmed].isSet("charge")) {
        const BoutReal charge = component_options[name_trimmed]["charge"];
        if (charge > 1e-5) {
          positive_ions.push_back(name_trimmed);
        } else if (charge < -1e-5) {
          negative_ions.push_back(name_trimmed);
        } else {
          neutrals.push_back(name_trimmed);
        }
      } else {
        neutrals.push_back(name_trimmed);
      }
    }

    // For each component e.g. "e", several Component types can be created
    // but if types are not specified then the component name is used
    const std::string types =
        component_options[name_trimmed].isSet("type")
            ? component_options[name_trimmed]["type"].as<std::string>()
            : name_trimmed;

    for (const auto &type : strsplit(types, ',')) {
      auto type_trimmed = trim(type, " \t\r()");
      if (type_trimmed.empty()) {
        continue;
      }

      components.push_back(Component::create(type_trimmed,
                                             name_trimmed,
                                             component_options,
                                             solver));
    }
  }

  const SpeciesInformation species(electrons, neutrals, positive_ions, negative_ions);
    
  for (auto& component : components) {
    component->declareAllSpecies(species);
  }

  components = sortComponents(std::move(components));
}

std::unique_ptr<ComponentScheduler> ComponentScheduler::create(Options &scheduler_options,
                                                               Options &component_options,
                                                               Solver *solver) {
  return std::make_unique<ComponentScheduler>(scheduler_options,
                                              component_options, solver);
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

void ComponentScheduler::outputVars(Options &state) {
  // Run through each component
  for(auto &component : components) {
    component->outputVars(state);
  }
}

void ComponentScheduler::restartVars(Options &state) {
  // Run through each component
  for(auto &component : components) {
    component->restartVars(state);
  }
}

void ComponentScheduler::precon(const Options &state, BoutReal gamma) {
  for(auto &component : components) {
    component->precon(state, gamma);
  }
}
