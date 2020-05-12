
#include "../include/component.hxx"

std::unique_ptr<Component> Component::create(const std::string &type,
                                             const std::string &name,
                                             Options &alloptions,
                                             Solver *solver) {

  return std::unique_ptr<Component>(
      ComponentFactory::getInstance().create(type, name, alloptions, solver));
}
