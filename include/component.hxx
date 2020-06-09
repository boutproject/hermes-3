#pragma once

#ifndef HERMES_COMPONENT_H
#define HERMES_COMPONENT_H

#include <options.hxx>
#include <bout/generic_factory.hxx>

#include <map>
#include <string>
#include <memory>

class Solver; // Time integrator

/// Interface for a component of a simulation model
/// 
/// The constructor of derived types should have signature
///   (std::string name, Options &options, Solver *solver)
/// 
struct Component {
  /// Modify the given simulation state
  virtual void transform(Options &state) = 0;
  
  /// Use the final simulation state to update internal state
  /// (e.g. time derivatives)
  virtual void finally(const Options &UNUSED(state)) { }

  /// Add extra fields for output, or set attributes e.g docstrings
  virtual void annotate(Options &UNUSED(state)) { }

  /// Preconditioning
  virtual void precon(const Options &UNUSED(state), BoutReal UNUSED(gamma)) { }
  
  /// Create a Component
  static std::unique_ptr<Component> create(const std::string &type, // The type to create
                                           const std::string &name, // The species/name for this instance
                                           Options &options,  // Component settings: options[name] are specific to this component
                                           Solver *solver); // Time integration solver
};

///////////////////////////////////////////////////////////////////

using ComponentCreator =
  std::function<Component *(std::string, Options &, Solver *)>;

class ComponentFactory : public Factory<Component, ComponentFactory, ComponentCreator> {
public:
  static constexpr auto type_name = "Component";
  static constexpr auto section_name = "component";
  static constexpr auto option_name = "type";
  static constexpr auto default_type = "none";
};

template <typename DerivedType>
class RegisterInFactory<Component, DerivedType, ComponentFactory> {
public:
  RegisterInFactory(const std::string &type) {
    ComponentFactory::getInstance().add(
        type,
        [](const std::string &name, Options &options, Solver *solver)
            -> Component * { return new DerivedType(name, options, solver); });
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
using RegisterComponent = RegisterInFactory<Component, DerivedType, ComponentFactory>;


/// Faster non-printing getter for Options
/// If this fails, it will throw BoutException
template<typename T>
T get(const Options& option) {
  if (!option.isSet()) {
    throw BoutException("Option {:s} has no value", option.str());
  }
  try {
    return bout::utils::variantStaticCastOrThrow<Options::ValueType, T>(option.value);
  } catch (const std::bad_cast &e) {
    // Convert to a more useful error message
    throw BoutException("Could not convert {:s} to type {:s}",
                        option.str(), typeid(T).name());
  }
}

/// Set values in an option. This could be optimised, but
/// currently the is_value private variable would need to be modified.
template<typename T>
Options& set(Options& option, T value) {
  option.force(std::move(value));
  return option;
}

/// Add value to a given option. If not already set, treats
/// as zero and sets the option to the value.
template<typename T>
Options& add(Options& option, T value) {
  if (!option.isSet()) {
    return set(option, value);
  } else {
    try {
      return set(option, value + bout::utils::variantStaticCastOrThrow<Options::ValueType, T>(option.value));
    } catch (const std::bad_cast &e) {
      // Convert to a more useful error message
      throw BoutException("Could not convert {:s} to type {:s}",
                          option.str(), typeid(T).name());
    }
  }
}

/// Add value to a given option. If not already set, treats
/// as zero and sets the option to the value.
template<typename T>
Options& subtract(Options& option, T value) {
  if (!option.isSet()) {
    option = -value;
  } else {
    try {
      set(option, bout::utils::variantStaticCastOrThrow<Options::ValueType, T>(option.value) - value);
    } catch (const std::bad_cast &e) {
      // Convert to a more useful error message
      throw BoutException("Could not convert {:s} to type {:s}",
                          option.str(), typeid(T).name());
    }
  }
  return option;
}

#endif // HERMES_COMPONENT_H
