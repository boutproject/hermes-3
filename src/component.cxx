
#include "../include/component.hxx"

Datafile *restart_datafile; ///< Temporary hack, to allow save/load from restarts

Datafile *get_restart_datafile() {
  return restart_datafile;
}

void set_restart_datafile(Datafile *file) {
  restart_datafile = file;
}

std::unique_ptr<Component> Component::create(const std::string &type,
                                             const std::string &name,
                                             Options &alloptions,
                                             Solver *solver) {

  return std::unique_ptr<Component>(
      ComponentFactory::getInstance().create(type, name, alloptions, solver));
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

bool isSetFinalNoBoundary(const Options& option, const std::string& location) {
#if CHECKLEVEL >= 1
  // Mark option as final inside the domain, but not in the boundary
  const_cast<Options&>(option).attributes["final-domain"] = location;
#endif
  return option.isSet();
}
