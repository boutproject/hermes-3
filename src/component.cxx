#include <memory>
#include <string>

#include <bout/assert.hxx>
#include <bout/options.hxx>
#include <bout/output.hxx>

#include "../include/component.hxx"
#include "../include/guarded_options.hxx"
#include "../include/permissions.hxx"

std::unique_ptr<Component> Component::create(const std::string &type,
                                             const std::string &name,
                                             Options &alloptions,
                                             Solver *solver) {

  return ComponentFactory::getInstance().create(type, name, alloptions, solver);
}

void Component::transform(Options& state) {
  GuardedOptions guarded(&state, &state_variable_access);
  transform_impl(guarded);
#if CHECKLEVEL >= 1
  for (auto& [varname, region] : guarded.unreadItems()) {
    output_warn.write("Did not read from state variable {} in region(s) {}\n", varname,
                      Permissions::regionNames(region));
  }
  for (auto& [varname, region] : guarded.unwrittenItems()) {
    output_warn.write("Did not write to state variable {} in region(s) {}\n", varname,
                      Permissions::regionNames(region));
  }
#endif
}

void Component::performAllSubstitutions(const SpeciesInformation& info) {
  // FIXME: doing all substitutions at once means that if any of these
  // species substitutions are empty then no substitutions are
  // performed.
  registerSubstitution("electrons", info.electrons);
  registerSubstitution("electrons2", info.electrons);
  registerSubstitution("neutrals", info.neutrals);
  registerSubstitution("neutrals2", info.neutrals);
  registerSubstitution("positive_ions", info.positive_ions);
  registerSubstitution("positive_ions2", info.positive_ions);
  registerSubstitution("negative_ions", info.negative_ions);
  registerSubstitution("negative_ions2", info.negative_ions);
  registerSubstitution("ions", info.ions);
  registerSubstitution("ions2", info.ions);
  registerSubstitution("charged", info.charged);
  registerSubstitution("charged2", info.charged);
  registerSubstitution("non_electrons", info.non_electrons);
  registerSubstitution("non_electrons2", info.non_electrons);
  registerSubstitution("all_species", info.all_species);
  registerSubstitution("all_species2", info.all_species);
  state_variable_access.substitute(permission_substitutions);
}

constexpr decltype(ComponentFactory::type_name) ComponentFactory::type_name;
constexpr decltype(ComponentFactory::section_name) ComponentFactory::section_name;
constexpr decltype(ComponentFactory::option_name) ComponentFactory::option_name;
constexpr decltype(ComponentFactory::default_type) ComponentFactory::default_type;

bool isSetFinal(const Options& option, const std::string& location) {
#if CHECKLEVEL >= 1
  // Mark option as final, both inside the domain and the boundary
  const_cast<Options&>(option).attributes["final"] = location;
  const_cast<Options&>(option).attributes["final-domain"] = location;
#endif
  return option.isSet();
}

bool isSetFinal(const GuardedOptions & option, const std::string& location) {
  const bool set = option.isSet();
#if CHECKLEVEL >= 1
  const PermissionTypes perm = option.getHighestPermission();
  if (perm >= PermissionTypes::Read
      or (perm == PermissionTypes::ReadIfSet and set)) {
    const Options& opt = option.get();
    const_cast<Options&>(opt).attributes["final"] = location;
    const_cast<Options&>(opt).attributes["final-domain"] = location;
  }
#endif
  return set;
}


bool isSetFinalNoBoundary(const Options& option, const std::string& location) {
#if CHECKLEVEL >= 1
  // Mark option as final inside the domain, but not in the boundary
  const_cast<Options&>(option).attributes["final-domain"] = location;
#endif
  return option.isSet();
}

bool isSetFinalNoBoundary(const GuardedOptions & option, const std::string& location) {
  const bool set = option.isSet();
#if CHECKLEVEL >= 1
  const PermissionTypes perm = option.getHighestPermission(Regions::Interior);
  if (perm >= PermissionTypes::Read or (perm == PermissionTypes::ReadIfSet and set)) {
    // Mark option as final inside the domain, but not in the boundary
    const_cast<Options&>(option.get(Regions::Interior)).attributes["final-domain"] =
        location;
  }
#endif
  return set;
}
