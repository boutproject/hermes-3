
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
