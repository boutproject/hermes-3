
#include "../include/component.hxx"

std::unique_ptr<Component> Component::create(const std::string &name,
                                             Options &options,
                                             const MeshMap &meshes) {
  return std::unique_ptr<Component>(ComponentFactory::getInstance().create(
      options["type"], name, options, meshes));
}
