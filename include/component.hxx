#pragma once

#ifndef HERMES_COMPONENT_H
#define HERMES_COMPONENT_H

#include <options.hxx>
#include <bout/generic_factory.hxx>

#include <map>
#include <string>
#include <memory>

class Mesh; // Probably already included in options

using MeshMap = std::map<std::string, Mesh*>;

/// Interface for a component of a simulation model
/// 
/// The constructor of derived types should have signature
///   (std::string name, Options &options, std::map<std::string, Mesh*> meshes)
/// 
struct Component {
  /// Modify the given simulation state
  virtual void transform(Options &state) = 0;
  
  /// Use the final simulation state to update internal state
  /// (e.g. time derivatives)
  virtual void finally(const Options &UNUSED(state)) { }
  
  /// Create a Component
  static std::unique_ptr<Component> create(const std::string &name, // The species/name for this instance
                                           Options &options,  // Settings from input
                                           const MeshMap &meshes); 
};

///////////////////////////////////////////////////////////////////

using ComponentCreator =
    std::function<Component *(std::string, Options &, const MeshMap &)>;

using ComponentFactory = Factory<Component, ComponentCreator>;

template <typename DerivedType>
class RegisterInFactory<Component, DerivedType> {
public:
  RegisterInFactory(const std::string &type) {
    ComponentFactory::getInstance().add(
        type,
        [](const std::string &name, Options &options, const MeshMap &meshes)
            -> Component * { return new DerivedType(name, options, meshes); });
  }
};

/// Simpler name for Factory registration helper class
///
/// Usage:
///
///     #include "component.hxx"
///     namespace {
///     RegisterComponent<MyComponent> registercomponentmine("mycomponent");
///     }
template <typename DerivedType>
using RegisterComponent = RegisterInFactory<Component, DerivedType>;


/// Faster non-printing getter for Options
/// If this fails, it will throw BoutException
template<typename T>
T get(const Options& option) {
  try {
    return bout::utils::variantStaticCastOrThrow<Options::ValueType, T>(option.value);
  } catch (const std::bad_cast &e) {
    // Convert to a more useful error message
    throw BoutException("Could not convert %s to type %s",
                        option.str().c_str(), typeid(T).name());
  }
}

/// Set values in an option. This could be optimised, but
/// currently the is_value private variable would need to be modified.
template<typename T>
Options& set(Options& option, T value) {
  option.force(std::move(value));
  return option;
}

#endif // HERMES_COMPONENT_H
