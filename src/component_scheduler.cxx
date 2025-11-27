#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <bout/bout_types.hxx>
#include <bout/options.hxx>
#include <bout/utils.hxx> // for trim, strsplit
#include <fmt/ranges.h>

#include "../include/component.hxx"
#include "../include/component_scheduler.hxx"

const std::set<std::string> ComponentScheduler::predeclared_variables = {
    "time",          "linear",      "units:inv_meters_cubed", "units:eV", "units:Tesla",
    "units:seconds", "units:meters"};

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

std::set<std::string> getParents(const std::string& name) {
  std::set<std::string> result;
  size_t start = 0, position = name.find(":", start);
  while (position != std::string::npos) {
    result.insert(name.substr(0, position));
    start = position + 1;
    position = name.find(":", start);
  }
  return result;
}

/// Produce a map between Option paths and all variable names held within
/// that path. If a path does not refer to a section then it just maps
/// to itself.
std::map<std::string, std::set<std::string>>
getVariableHierarchy(const std::vector<std::unique_ptr<Component>>& components) {
  // Build up a set of all variable names which are read only if they
  // are set by another component
  std::set<std::string> conditional_names;
  for (const auto& component : components) {
    const Permissions& permissions = component->getPermissions();
    for (const auto& [varname, _] :
         permissions.getVariablesWithPermission(PermissionTypes::ReadIfSet, true)) {
      conditional_names.insert(varname);
    }
  }

  // Build up a set of all variable names which are definitely
  // read/written by components, and the sections which they imply
  // exist
  std::set<std::string> unconditional_names, unconditional_sections;
  for (const auto& component : components) {
    const Permissions& permissions = component->getPermissions();
    for (const auto& [varname, _] :
         permissions.getVariablesWithPermission(PermissionTypes::Read, false)) {
      unconditional_names.insert(varname);
      unconditional_sections.merge(getParents(varname));
    }
  }

  // Split the set of all variable names used by components into
  // those that are sections and those that are not.
  std::set<std::string> sections_present, non_sections;
  std::set_intersection(unconditional_names.begin(), unconditional_names.end(),
                        unconditional_sections.begin(), unconditional_sections.end(),
                        std::inserter(sections_present, sections_present.begin()));
  std::set_difference(unconditional_names.begin(), unconditional_names.end(),
                      sections_present.begin(), sections_present.end(),
                      std::inserter(non_sections, non_sections.begin()));

  std::map<std::string, std::set<std::string>> result;

  // ReadIfSet variables will only actually be used if they are
  // reference elsewhere. We create them with empty sets, which will
  // get filled if they are present.
  for (const auto& name : conditional_names) {
    result.insert({name, {name}});
  }

  // Non-sections map to themselves
  for (const auto& name : non_sections) {
    result[name] = {name};
  }
  // Sections map to those variables which they contain
  for (const auto& section : sections_present) {
    auto& children = result[section];
    const std::string sec_suffixed = section + ':';
    for (const auto& name : non_sections) {
      if (name.rfind(sec_suffixed, 0) == 0) {
        children.insert(name);
      }
    }
  }

  return result;
}

/// Get all variables to which the name could be referring (e.g., its
/// children if it is a section name). These will be filtered to
/// remove any variables for which more specific permissions are
/// given.
std::set<std::string>
expandVariableName(const std::map<std::string, std::set<std::string>>& hierarchy,
                   const Permissions& permissions, const std::string& name) {
  const std::set<std::string>& candidates = hierarchy.at(name);
  std::set<std::string> result;
  // Only return the values that do not have a more specific permission
  std::copy_if(candidates.begin(), candidates.end(),
               std::inserter(result, result.begin()),
               [&permissions, &name](const std::string& candidate) -> bool {
                 return permissions.bestMatchRights(candidate).name == name;
               });
  return result;
}

using Var = std::pair<std::string, Regions>;

std::map<Var, std::set<size_t>>
getPermissionComponentMap(const std::vector<std::unique_ptr<Component>>& components,
                          const std::map<std::string, std::set<std::string>>& hierarchy,
                          PermissionTypes permission) {
  std::map<Var, std::set<size_t>> result;
  for (size_t i = 0; i < components.size(); i++) {
    const Permissions& permissions = components[i]->getPermissions();
    for (const auto& [name, regions] :
         permissions.getVariablesWithPermission(permission)) {
      for (const auto& sub_name : expandVariableName(hierarchy, permissions, name)) {
        for (const auto& [region, _] : Permissions::fundamental_regions) {
          if ((regions & region) == region) {
            result[{sub_name, region}].insert(i);
          }
        }
      }
    }
  }
  return result;
}

/// Modifies component_dependencies to include information on which
/// variables each component reads. Returns a set of any read
/// variables which are not written by any component.
std::set<std::string>
setReadDependencies(const std::vector<std::unique_ptr<Component>>& components,
                    const std::map<std::string, std::set<std::string>>& hierarchy,
                    const std::map<Var, std::set<size_t>>& writers,
                    PermissionTypes permission,
                    std::vector<std::set<size_t>>& component_dependencies) {
  std::set<std::string> missing;
  for (size_t i = 0; i < components.size(); i++) {
    const Permissions& permissions = components[i]->getPermissions();
    // Create dependencies between components that read variables and those that write
    // them
    for (const auto& [name, regions] :
         permissions.getVariablesWithPermission(permission)) {
      if (ComponentScheduler::predeclared_variables.count(name) > 0)
        continue;
      for (const auto& sub_name : expandVariableName(hierarchy, permissions, name)) {
        for (const auto& [region, _] : Permissions::fundamental_regions) {
          if ((regions & region) == region) {
            const auto item = writers.find({sub_name, region});
            if (item == writers.end()) {
              missing.insert(
                  fmt::format("{} ({})", sub_name, Permissions::regionNames(region)));
            } else {
              component_dependencies[i].insert(item->second.begin(), item->second.end());
            }
          }
        }
      }
    }
  }
  return missing;
}

/// Consumes a list of components and returns a new one that has been
/// topolgically sorted to ensure variables are written and read in
/// the right order.
std::vector<std::unique_ptr<Component>>
sortComponents(std::vector<std::unique_ptr<Component>>&& components) {
  // Map variables to the components that write them
  std::map<std::string, std::set<std::string>> variable_hierarchy =
      getVariableHierarchy(components);

  // Get information on which components write each variable
  std::map<Var, std::set<size_t>> nonfinal_writes = getPermissionComponentMap(
                                      components, variable_hierarchy,
                                      PermissionTypes::Write),
                                  final_writes = getPermissionComponentMap(
                                      components, variable_hierarchy,
                                      PermissionTypes::Final);

  std::vector<std::set<size_t>> component_dependencies(components.size());

  // Components which do a final write on a variable depend on all
  // components which do non-final writes on that variable
  for (const auto& [var, comp_indices] : final_writes) {
    for (size_t i : comp_indices) {
      const auto item = nonfinal_writes.find(var);
      if (item != nonfinal_writes.end()) {
        // Note that calling merge actually removes the items from
        // nonfinal_writes. This is fine because the only remaining
        // thing for which we will use nonfinal_writes is setting up
        // variable_writers. That doesn't use any information on
        // nonfinal writes for variables which have a final write, so
        // it won't do any harm to remove it.
        component_dependencies[i].merge(item->second);
      }
    }
  }

  // Work out which component(s) last write a variable before it may be read
  std::map<Var, std::set<size_t>> variable_writers = std::move(final_writes);
  variable_writers.merge(std::move(nonfinal_writes));

  // Create dependency information for read variables
  std::set<std::string> missing =
      setReadDependencies(components, variable_hierarchy, variable_writers,
                          PermissionTypes::Read, component_dependencies);
  if (missing.size() > 0) {
    throw BoutException(
        "The following required variables are not written by any component:\n\t{}\n",
        fmt::format("{}", fmt::join(missing, "\n\t")));
  }
  // If can not find a place where a ReadIfSet variable is written, it will just be
  // skipped
  setReadDependencies(components, variable_hierarchy, variable_writers,
                      PermissionTypes::ReadIfSet, component_dependencies);

  // Create ancillary variables for sorting process
  std::vector<bool> processing(components.size(), false),
      processed(components.size(), false);
  std::vector<std::size_t> order;

  // Perform the sort
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
