#include <memory>
#include <string>
#include <vector>

#include <bout/bout_types.hxx>
#include <bout/options.hxx>
#include <bout/utils.hxx> // for trim, strsplit

#include "../include/component.hxx"
#include "../include/component_scheduler.hxx"

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
}

// Use index to identify component
// std::map<int, std::set<int>> dependencies;
// std::map<std::string, int> var_written_by, var_final_written_by, read_dependencies;

// Assemble the maps of which components write each variable
// Assemble read_dependencies using values from var_final_written_by, if present, else
// from var_written_by For each component
//   For each write-final var, create dependency between component and all components that
//   write it For each read-if-set var that has a write-final, create dependency between
//   component and (final) writer(s) For each read-var create a dependency between
//   component and (final) writers; if none then raise exception

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
