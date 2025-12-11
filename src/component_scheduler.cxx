#include <algorithm>
#include <cstddef>
#include <iterator>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/options.hxx>
#include <bout/utils.hxx> // for trim, strsplit
#include <fmt/format.h>
#include <fmt/ranges.h>

#include "../include/component.hxx"
#include "../include/component_scheduler.hxx"
#include "../include/permissions.hxx"

const std::set<std::string> ComponentScheduler::predeclared_variables = {
    "time",          "linear",      "units:inv_meters_cubed", "units:eV", "units:Tesla",
    "units:seconds", "units:meters"};

/// Perform a depth-first topological sort, starting from `item`. It
/// will finish once it reaches the end of `item`'s dependency chain,
/// so this needs to be called in a loop for all items. Information
/// about `item` is stored in the corresponding index of the vector
/// arguments.
///
/// In pratice, `item` here represents the index of a particular
/// component. The indices of the dependencies of `item` are stored in
/// the corresponding element of `dependencies`.
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

/// Get all the parent sections of a variable "path". Sections are
/// separated by colons in the path.
std::set<std::string> getParents(const std::string& name) {
  std::set<std::string> result;
  size_t start = 0;
  size_t position = name.find(":", start);
  while (position != std::string::npos) {
    result.insert(name.substr(0, position));
    start = position + 1;
    position = name.find(":", start);
  }
  return result;
}

/// Produce a map between Option paths and all variable names held
/// within that path. If the path refers to a section then it maps to
/// the set of all variables contained in that section and any
/// sub-sections. Otherwise the path corresponds to a variable and
/// just maps to itself. Only paths which are explicitly given a
/// permission by at least one component will be present.
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

  // Build up a set of all section/variable names which are definitely
  // read/written by components, and the sections which they imply
  // exist
  std::set<std::string> unconditional_names;
  std::set<std::string> unconditional_sections;
  for (const auto& component : components) {
    const Permissions& permissions = component->getPermissions();
    for (const auto& [varname, _] :
         permissions.getVariablesWithPermission(PermissionTypes::Read, false)) {
      unconditional_names.insert(varname);
      unconditional_sections.merge(getParents(varname));
    }
  }

  /// Assemble the set of all section names which are referred to
  /// explicitly in the component permissions.
  std::set<std::string> sections_present;
  std::set_intersection(unconditional_names.begin(), unconditional_names.end(),
                        unconditional_sections.begin(), unconditional_sections.end(),
                        std::inserter(sections_present, sections_present.begin()));
  /// Assemble the set of all variable names which are definitlye
  /// read/written by components (i.e., not including sections)
  std::set<std::string> non_sections;
  std::set_difference(unconditional_names.begin(), unconditional_names.end(),
                      sections_present.begin(), sections_present.end(),
                      std::inserter(non_sections, non_sections.begin()));

  std::map<std::string, std::set<std::string>> result;

  // ReadIfSet variables will only actually be used if they are
  // reference elsewhere. We create them with empty sets, which will
  // get filled if they are present.
  for (const auto& name : conditional_names) {
    // FIXME: this isn't an empty set
    //result.insert({name, {name}});
    result.insert({name, {}});
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

/// Create a map between a variable and the set of components that
/// access it with the specified permission level.
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
/// components depend on each other. It does this by making components
/// which read a variable depend on whichever component(s) write that
/// variable (information contained in the `writers`
/// argument). Returns a set of any read variables which are not
/// written by any component.
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
      if (ComponentScheduler::predeclared_variables.count(name) > 0) {
        continue;
      }
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

/// Topologically sorts the list of components to ensure variables are
/// written and read in the right order.
void sortComponents(std::vector<std::unique_ptr<Component>>& components) {
  // Map between variable/section names specified by component
  // permissions and the variables they contain. In the case of
  // sections this is all variables within the section and any
  // sub-sections. Non-section viarables map to themselves.
  const std::map<std::string, std::set<std::string>> variable_hierarchy =
      getVariableHierarchy(components);

  // Get information on which components write each variable
  std::map<Var, std::set<size_t>> nonfinal_writes =
      getPermissionComponentMap(components, variable_hierarchy, PermissionTypes::Write);
  std::map<Var, std::set<size_t>> final_writes =
      getPermissionComponentMap(components, variable_hierarchy, PermissionTypes::Final);

  // Object mapping between components (reprsented by the index of
  // that component in the `components` argument) and the components
  // each of these depends upon (represented by a set of the indices
  // for those components).
  std::vector<std::set<size_t>> component_dependencies(components.size());

  // Components which do a final write on a variable depend on all
  // components which do non-final writes on that variable
  for (const auto& [var, comp_indices] : final_writes) {
    if (comp_indices.size() > 1) {
      throw BoutException(
          "Multiple components have permission to make final write to variable {}", var);
    }
    for (const size_t i : comp_indices) {
      const auto item = nonfinal_writes.find(var);
      if (item != nonfinal_writes.end()) {
        // Note that calling merge actually removes the items from the
        // sets stored in nonfinal_writes. This is fine because the
        // only remaining thing for which we will use nonfinal_writes
        // is setting up variable_writers and that doesn't use any
        // information on variables which have a final
        // write. Therefore it won't do any harm hear to remove
        // information about variables which have a final write..
        component_dependencies[i].merge(item->second);
      }
    }
  }

  // Work out which component(s) last write a variable before it may
  // be read. For variables with a final-write, it is whichever
  // component performs that final write.
  std::map<Var, std::set<size_t>> variable_writers = std::move(final_writes);
  // For other variables, it is the set of all components which have
  // write permission.
  variable_writers.merge(std::move(nonfinal_writes));

  // Insert dependency information for components that (unconditionally) read variables
  std::set<std::string> missing =
      setReadDependencies(components, variable_hierarchy, variable_writers,
                          PermissionTypes::Read, component_dependencies);
  if (missing.size() > 0) {
    throw BoutException(
        "The following required variables are not written by any component:\n\t{}\n",
        fmt::format("{}", fmt::join(missing, "\n\t")));
  }
  // Insert dependency information for components that read variables
  // if those variables have been set. If can not find a place where
  // the variable is written, it will just be skipped.
  setReadDependencies(components, variable_hierarchy, variable_writers,
                      PermissionTypes::ReadIfSet, component_dependencies);

  // Create ancillary variables for sorting process
  std::vector<bool> processing(components.size(), false);
  std::vector<bool> processed(components.size(), false);
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

  components = std::move(result);
}

ComponentScheduler::ComponentScheduler(Options &scheduler_options,
                                       Options &component_options,
                                       Solver *solver) {

  const std::string component_names = scheduler_options["components"]
                                          .doc("Components in order of execution")
                                          .as<std::string>();

  std::set<ComponentInformation> required_components;

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

      required_components.emplace(name_trimmed, type_trimmed);
    }
  }

  std::set<ComponentInformation> created_components;
  for (auto it = required_components.begin(); it != required_components.end();
       it = required_components.begin()) {
    auto comp = Component::create(it->type, it->name, component_options, solver);
    for (const auto& sub_comp : comp->additionalComponents()) {
      if (required_components.count(sub_comp) == 0
          and created_components.count(sub_comp) == 0) {
        required_components.insert(sub_comp);
      }
    }
    components.push_back(std::move(comp));
    created_components.insert(required_components.extract(it));
  }

  const SpeciesInformation species(electrons, neutrals, positive_ions, negative_ions);
    
  for (auto& component : components) {
    component->declareAllSpecies(species);
  }

  sortComponents(components);
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
