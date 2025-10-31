
#include "../include/component.hxx"

std::unique_ptr<Component> Component::create(const std::string &type,
                                             const std::string &name,
                                             Options &alloptions,
                                             Solver *solver) {

  return ComponentFactory::getInstance().create(type, name, alloptions, solver);
}

void Component::transform(Options& state) {
  GuardedOptions guarded(&state, &state_variable_access);
  transform(guarded);
  for (auto& [varname, region] : guarded.unreadItems()) {
    output_warn.write("Did not read from state variable {} in region(s) {}",
                      varname, Permissions::regionNames(region));
  }
  for (auto& [varname, region] : guarded.unwrittenItems()) {
    output_warn.write("Did not write to state variable {} in region(s) {}",
                      varname, Permissions::regionNames(region));
  }
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

bool isSetFinal(const GuardedOptions option, const std::string& location) {
  return isSetFinal(option.get(), location);
}


bool isSetFinalNoBoundary(const Options& option, const std::string& location) {
#if CHECKLEVEL >= 1
  // Mark option as final inside the domain, but not in the boundary
  const_cast<Options&>(option).attributes["final-domain"] = location;
#endif
  return option.isSet();
}

bool isSetFinalNoBoundary(const GuardedOptions option, const std::string& location) {
  return isSetFinalNoBoundary(option.get(Permissions::Interior), location);
}
