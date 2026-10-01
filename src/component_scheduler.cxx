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
#include <bout/output.hxx>
#include <bout/utils.hxx> // for trim, strsplit
#include <fmt/format.h>
#include <fmt/ranges.h>

#include "../include/component.hxx"
#include "../include/component_scheduler.hxx"
#include "../include/permissions.hxx"

namespace {

/// Cheap class to track dependencies while recursing through the
/// dependency graph. It is designed to minimise the ammount of
/// copying that needs to be done.
class DependencyChain {
private:
  struct Element {
    Element(std::shared_ptr<Element> _next, Component* _component)
        : next(std::move(_next)), component(_component) {}
    std::shared_ptr<Element> next;
    Component* component;
  };

public:
  DependencyChain(std::vector<std::unique_ptr<Component>>* components)
      : components(components) {}

  DependencyChain addComponent(Component* component) const {
    return DependencyChain(*this, component);
  }
  DependencyChain addComponent(const size_t comp_idx) const {
    return addComponent((*components)[comp_idx].get());
  }

  struct iterator {
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = Component;
    using pointer = Component*;
    using reference = Component&;

    iterator(std::shared_ptr<Element> item) : current(std::move(item)) {}

    reference operator*() const { return *current->component; }
    pointer operator->() { return current->component; }

    iterator& operator++() {
      current = current->next;
      return *this;
    }
    iterator operator++(int) {
      iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    friend bool operator==(const iterator& a, const iterator& b) {
      return a.current == b.current;
    }
    friend bool operator!=(const iterator& a, const iterator& b) {
      return a.current != b.current;
    }

  private:
    std::shared_ptr<Element> current;
  };

  iterator begin() const { return back; }
  iterator end() const { return std::shared_ptr<Element>(nullptr); }

  std::size_t size() const { return _size; }

private:
  std::vector<std::unique_ptr<Component>>* components = nullptr;
  std::shared_ptr<Element> back{nullptr};
  std::size_t _size = 0;

  DependencyChain(const DependencyChain& tail, Component* new_comp)
      : components(tail.components), back(std::make_shared<Element>(tail.back, new_comp)),
        _size(tail._size + 1) {}
};

/// All the variable names which are pre-set in the state, before
/// any components are applied.
const std::set<std::string> predeclared_variables = {
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
                      std::vector<bool>& processed, const DependencyChain& dependees) {
  if (processed[item]) {
    return;
  }
  if (processing[item]) {
    auto recursed_dependees = dependees.addComponent(item);
    throw BoutException("Circular dependency among components: {}",
                        fmt::join(recursed_dependees, " -> "));
  }
  processing[item] = true;

  if (not dependencies[item].empty()) {
    auto recursed_dependees = dependees.addComponent(item);
    for (const auto dep : dependencies[item]) {
      topological_sort(dependencies, dep, sorted, processing, processed,
                       recursed_dependees);
    }
  }
  processed[item] = true;
  sorted.push_back(item);
}

/// Get all the parent sections of a variable "path" (i.e., in the
/// hierarchy of `Options` objects in the state). Sections are
/// separated by colons in the path. This is different from just
/// splitting on ':' because it returns the hierarchy of
/// fully-qualified parent sections. E.g.
///
///     getParents("species:d:collision_frequencies:d_d_coll")
///
/// would return {"species", "species:d",
/// "species:d:collision_frequencies"}.
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
/// the set of all variables contained in that section and its
/// sub-sections. Otherwise the path corresponds to a variable and
/// just maps to itself. Only paths which are explicitly given a
/// permission by at least one component will be present. Paths with
/// only `readIfSet` permission will only map to anything if there is
/// a component that has write permission for it or one of its
/// parents.
///
/// The algorithm for this is:
///
///  - Construct a set of all names that have the following permissions. These names may refer
///    either to variable names or section names; at this stage we don't know which.
///    - readIfSet permission (readifset_names)
///    - read, write, or writeFinal permission (readwrite_names)
///    - write or writeFinal permission (write_names)
///
///  - Assemble sets consisting of the names of the parent sections
///    for the contents of readifset_names (readifset_sections) and
///    readwrite_names (readwrite_sections)
///
///  - Use the above to build the following further sets:
///    - All variable/section names with permissions
///        all_names = readifset_names ∪ readwrite_names
///    - The parent sections of the contents of all_names
///        all_sections = readifset_sections ∪ readwrite_sections
///    - Section names explicitly given readIfSet permissions
///        readifset_sections_present = readifset_names ∩ all_sections
///    - Variable names explicitly given readIfSet permissions
///        readifset_non_sections = readifset_names \ readifset_sections_present
///    - Section names explicitly given read permissions or higher
///        readwrite_sections_present = readwrite_names ∩ all_sections
///    - Variable names explicitly given read permissions or higher
///        readwrite_non_sections = readwrite_names \ readwrite_sections_present
///    - Section names explicitly given write permissions or higher
///        write_sections_present = write_names ∩ all_sections
///    - Variable names explicitly given write permissions or higher
///         write_non_sections = write_names \ write_sections_present
///
///  - Build up the map between names and variables using the following rules:
///    - The contents of readifset_non_sections map to themselves if that variable
///      name or one of its parents has been given write permission elsewhere
///    - The contents of readifset_sections_present map to any children variables
///      given write permission elsewhere
///    - The contents of readwrite_non_sections map to themselves
///    - The contents of readwrite_sections_present map to any children variables
///      given read or write permission elsewhere.
///    - Additionally, the contents of write_sections_present will map to any child
///      variables which have readIfSet permissions.
std::map<std::string, std::set<std::string>>
getVariableHierarchy(const std::vector<std::unique_ptr<Component>>& components) {
  std::set<std::string> readifset_names;
  std::set<std::string> readifset_sections;
  std::set<std::string> readwrite_names;
  std::set<std::string> readwrite_sections;
  std::set<std::string> write_names;
  for (const auto& component : components) {
    const Permissions& permissions = component->getPermissions();
    // Build up a set of all variable names which are read only if they
    // are set by another component
    for (const auto& [varname, _] :
         permissions.getVariablesWithPermission(PermissionTypes::ReadIfSet)) {
      readifset_names.insert(varname);
      readifset_sections.merge(getParents(varname));
    }
    // Build up a set of all section/variable names which are definitely
    // read/written by components, and the sections which they imply
    // exist
    for (const auto& [varname, _] :
         permissions.getVariablesWithMinimumPermission(PermissionTypes::Read)) {
      readwrite_names.insert(varname);
      readwrite_sections.merge(getParents(varname));
    }
    // Build up a set of all section/variable names which are
    // written by components, and the sections which they imply
    // exist
    for (const auto& [varname, _] :
         permissions.getVariablesWithMinimumPermission(PermissionTypes::Write)) {
      write_names.insert(varname);
    }
  }

  /// Assemble the list of all section/variable names which are explicitly given permissions
  std::set<std::string> all_names;
  std::set_union(readifset_names.begin(), readifset_names.end(), readwrite_names.begin(),
                 readwrite_names.end(), std::inserter(all_names, all_names.begin()));
  /// Assemble the set of all sections which the names in all_names imply exist
  std::set<std::string> all_sections;
  std::set_union(readifset_sections.begin(), readifset_sections.end(),
                 readwrite_sections.begin(), readwrite_sections.end(),
                 std::inserter(all_sections, all_sections.begin()));

  /// Assemble the set of all section names which are explicitly given readIfSet permissions.
  std::set<std::string> readifset_sections_present;
  std::set_intersection(
      readifset_names.begin(), readifset_names.end(), all_sections.begin(),
      all_sections.end(),
      std::inserter(readifset_sections_present, readifset_sections_present.begin()));
  /// Assemble the set of all variable names which are given readIfSet
  /// permission and which are not sections
  std::set<std::string> readifset_non_sections;
  std::set_difference(
      readifset_names.begin(), readifset_names.end(), readifset_sections_present.begin(),
      readifset_sections_present.end(),
      std::inserter(readifset_non_sections, readifset_non_sections.begin()));

  /// Assemble the set of all section names which are explicitely
  /// given read permission or higher.
  std::set<std::string> readwrite_sections_present;
  std::set_intersection(
      readwrite_names.begin(), readwrite_names.end(), all_sections.begin(),
      all_sections.end(),
      std::inserter(readwrite_sections_present, readwrite_sections_present.begin()));
  /// Assemble the set of all variable names which are definitely
  /// read/written by components and which are not sections
  std::set<std::string> readwrite_non_sections;
  std::set_difference(
      readwrite_names.begin(), readwrite_names.end(), readwrite_sections_present.begin(),
      readwrite_sections_present.end(),
      std::inserter(readwrite_non_sections, readwrite_non_sections.begin()));

  /// Assemble the set of all section names which are explicitely
  /// given write permission or higher.
  std::set<std::string> write_sections_present;
  std::set_intersection(
      write_names.begin(), write_names.end(), all_sections.begin(), all_sections.end(),
      std::inserter(write_sections_present, write_sections_present.begin()));
  /// Assemble the set of all variable names which are given write
  /// permission or higher and which are not sections
  std::set<std::string> write_non_sections;
  std::set_difference(write_names.begin(), write_names.end(),
                      write_sections_present.begin(), write_sections_present.end(),
                      std::inserter(write_non_sections, write_non_sections.begin()));

  std::map<std::string, std::set<std::string>> result;

  // ReadIfSet variables will be used if they or a parent section have
  // write permission somewhere. Only map them to themselves if that
  // is the case. Ensure that parent sections with write permission
  // will also map to the readIfSet variable.
  for (const auto& name : readifset_non_sections) {
    auto& val = result[name];
    if (write_non_sections.contains(name)) {
      val.insert(name);
    }
    for (const auto& parent : getParents(name)) {
      if (write_sections_present.contains(parent)) {
        val.insert(name);
        result[parent].insert(name);
      }
    }
  }
  // ReadIfSet sections map to any children variables that have write permission
  for (const auto& section : readifset_sections_present) {
    auto& children = result[section];
    const std::string sec_suffixed = section + ':';
    for (const auto& name : write_non_sections) {
      if (name.rfind(sec_suffixed, 0) == 0) {
        children.insert(name);
      }
    }
  }
  // Non-sections map to themselves
  for (const auto& name : readwrite_non_sections) {
    result[name] = {name};
  }
  // Sections map to those variables which they contain
  for (const auto& section : readwrite_sections_present) {
    auto& children = result[section];
    const std::string sec_suffixed = section + ':';
    for (const auto& name : readwrite_non_sections) {
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
      if (predeclared_variables.count(name) > 0) {
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

void printComponents(const std::vector<std::unique_ptr<Component>>& components) {
  if (!components.empty()) {
    output_info << "\nComponents will be executed in the following order:\n";
  }
  for (const auto& comp : components) {
    output_info << fmt::format("\t{}\n", *comp);
  }
  if (!components.empty()) {
    output_info << "\n";
  }
}

/// Topologically sorts the list of components to ensure variables are
/// written and read in the right order.
///
/// This is quite a complicated process. The steps are:
///
/// 1. Construct a map between names and the variable(s) to which
///    they refer (section names refer to all variables contained within
///    the section). This is used so that, when a permission is set for
///    a whole section, we can work out what are the actual variables
///    the permission applies to.
/// 2. Identify the components which have permission to do final writes
///    and non-final writes on each variable.
/// 3. Construct a map between variable names and which components last
///    write to them.
///    - For variables where a component has final write permission,
///      it is that component.
///    - Otherwise, it is all components which have non-final write
///      permission for the variable.
/// 4. Establish the dependencies between components.
///    - Components which have final-write permission for a variable
///      depend on any components that have non-final write permission for
///      that variable.
///    - Components that have read permission for a variable will depend
///      on whichever component(s) last write to that variable, as
///      determined in the previous step.
/// 5. Use this dependency information to perform a topological sort on
///    the components.
void sortComponents(std::vector<std::unique_ptr<Component>>& components) {
  // Map between variable/section names specified by component
  // permissions and the variables they contain. In the case of
  // sections this is all variables within the section and its
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
      std::vector<std::string> comps;
      for (const auto i : comp_indices) {
        comps.push_back(fmt::format("{}", *components[i]));
      }
      throw BoutException(
          "Multiple components have permission to make final write to variable {}: {}",
          var, fmt::join(comps, ", "));
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
  if (!missing.empty()) {
    throw BoutException(
        "The following required variables are not written by any component:\n\t{}\n",
        fmt::join(missing, "\n\t"));
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
      topological_sort(component_dependencies, i, order, processing, processed,
                       DependencyChain(&components));
    }
  }

  // Create the result with components in the desired order
  std::vector<std::unique_ptr<Component>> result(components.size());
  for (size_t i = 0; i < components.size(); i++) {
    std::swap(result[i], components[order[i]]);
  }

  components = std::move(result);
}
} // namespace

ComponentScheduler::ComponentScheduler(Options& scheduler_options,
                                       Options& component_options, Solver* solver) {

  const std::string component_names = scheduler_options["components"]
                                          .doc("Components in order of execution")
                                          .as<std::string>();
  const bool autosort = scheduler_options["autosort"]
                            .doc("Perform a topological sort to ensure components "
                                 "executed in the right order?")
                            .withDefault<bool>(true);

  std::vector<std::string> electrons;
  std::vector<std::string> neutrals;
  std::vector<std::string> positive_ions;
  std::vector<std::string> negative_ions;

  // For now split on ','. Something like "->" might be better
  for (const auto& name : strsplit(component_names, ',')) {
    // Ignore brackets, to allow these to be used to span lines.
    // In future brackets may be useful for complex scheduling

    auto name_trimmed = trim(name, " \t\r()");
    if (name_trimmed.empty()) {
      continue;
    }

    if (name_trimmed == "e" or name_trimmed == "ebeam") {
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

    for (const auto& type : strsplit(types, ',')) {
      auto type_trimmed = trim(type, " \t\r()");
      if (type_trimmed.empty()) {
        continue;
      }

      components.push_back(
          Component::create(type_trimmed, name_trimmed, component_options, solver));
    }
  }

  const SpeciesInformation species(electrons, neutrals, positive_ions, negative_ions);

  for (auto& component : components) {
    component->declareAllSpecies(species);
  }

  if (autosort) {
    ::sortComponents(components);
  }
  printComponents(components);
}

std::unique_ptr<ComponentScheduler> ComponentScheduler::create(Options& scheduler_options,
                                                               Options& component_options,
                                                               Solver* solver) {
  return std::make_unique<ComponentScheduler>(scheduler_options, component_options,
                                              solver);
}

void ComponentScheduler::transform(Options& state) {
  // Run through each component
  for (auto& component : components) {
    component->transform(state);
  }
  // Enable components to update themselves based on the final state
  for (auto& component : components) {
    component->finally(state);
  }
}

void ComponentScheduler::outputVars(Options& state) {
  // Run through each component
  for (auto& component : components) {
    component->outputVars(state);
  }
}

void ComponentScheduler::restartVars(Options& state) {
  // Run through each component
  for (auto& component : components) {
    component->restartVars(state);
  }
}

void ComponentScheduler::precon(const Options& state, BoutReal gamma) {
  for (auto& component : components) {
    component->precon(state, gamma);
  }
}
