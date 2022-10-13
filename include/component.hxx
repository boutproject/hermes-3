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
  virtual ~Component() {}

  /// Modify the given simulation state
  /// All components must implement this function
  virtual void transform(Options &state) = 0;
  
  /// Use the final simulation state to update internal state
  /// (e.g. time derivatives)
  virtual void finally(const Options &UNUSED(state)) { }

  /// Add extra fields for output, or set attributes e.g docstrings
  virtual void outputVars(Options &UNUSED(state)) { }

  /// Preconditioning
  virtual void precon(const Options &UNUSED(state), BoutReal UNUSED(gamma)) { }
  
  /// Create a Component
  ///
  /// @param type     The name of the component type to create (e.g. "evolve_density")
  /// @param name     The species/name for this instance.
  /// @param options  Component settings: options[name] are specific to this component
  /// @param solver   Time-integration solver
  static std::unique_ptr<Component> create(const std::string &type, // The type to create
                                           const std::string &name, // The species/name for this instance
                                           Options &options,  // Component settings: options[name] are specific to this component
                                           Solver *solver); // Time integration solver
};

///////////////////////////////////////////////////////////////////

/// A factory for creating Components on demand, based on a string type name
/// The template arguments after ComponentFactory are the types of the arguments
/// to the Component constructor.
class ComponentFactory
    : public Factory<Component, ComponentFactory, const std::string&, Options&, Solver*> {
public:
  static constexpr auto type_name = "Component";
  static constexpr auto section_name = "component";
  static constexpr auto option_name = "type";
  static constexpr auto default_type = "none";
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
using RegisterComponent = ComponentFactory::RegisterInFactory<DerivedType>;

/// Faster non-printing getter for Options
/// If this fails, it will throw BoutException
///
/// This version allows the value to be modified later
/// i.e. the value returned is not the "final" value.
///
/// @tparam T  The type the option should be converted to
///
/// @param option  The Option whose value will be returned
template<typename T>
T getNonFinal(const Options& option) {
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

#define TOSTRING_(x) #x
#define TOSTRING(x) TOSTRING_(x)


/// Faster non-printing getter for Options
/// If this fails, it will throw BoutException
///
/// This marks the value as final, both in the domain and the boundary.
/// Subsequent calls to "set" this option will raise an exception.
///
/// @tparam T  The type the option should be converted to
///
/// @param option  The Option whose value will be returned
/// @param location  An optional string to indicate where this value is used
template<typename T>
T get(const Options& option, const std::string& location = "") {
#if CHECKLEVEL >= 1
  // Mark option as final, both inside the domain and the boundary
  const_cast<Options&>(option).attributes["final"] = location;
  const_cast<Options&>(option).attributes["final-domain"] = location;
#endif
  return getNonFinal<T>(option);
}

/// Check if an option can be fetched
/// Sets the final flag so setting the value
/// afterwards will lead to an error
bool isSetFinal(const Options& option, const std::string& location = "");

#if CHECKLEVEL >= 1
/// A wrapper around isSetFinal() which captures debugging information
///
/// Usage:
///   if (IS_SET(option["value"]));
#define IS_SET(option) \
  isSetFinal(option, __FILE__ ":" TOSTRING(__LINE__))
#else
#define IS_SET(option) \
  isSetFinal(option)
#endif

/// Check if an option can be fetched
/// Sets the final flag so setting the value in the domain
/// afterwards will lead to an error
bool isSetFinalNoBoundary(const Options& option, const std::string& location = "");

#if CHECKLEVEL >= 1
/// A wrapper around isSetFinalNoBoundary() which captures debugging information
///
/// Usage:
///   if (IS_SET_NOBOUNDARY(option["value"]));
#define IS_SET_NOBOUNDARY(option) \
  isSetFinalNoBoundary(option, __FILE__ ":" TOSTRING(__LINE__))
#else
#define IS_SET_NOBOUNDARY(option) \
  isSetFinalNoBoundary(option)
#endif

#if CHECKLEVEL >= 1
/// A wrapper around get<>() which captures debugging information
///
/// Usage:
///   auto var = GET_VALUE(Field3D, option["value"]);
#define GET_VALUE(Type, option) \
  get<Type>(option, __FILE__ ":" TOSTRING(__LINE__))
#else
#define GET_VALUE(Type, option) \
  get<Type>(option)
#endif

/// Faster non-printing getter for Options
/// If this fails, it will throw BoutException
///
/// This marks the value as final in the domain.
/// The caller is assuming that the boundary values are non-final or invalid.
/// Subsequent calls to "set" this option will raise an exception,
/// but calls to "setBoundary" will not.
///
/// @tparam T  The type the option should be converted to
///
/// @param option  The Option whose value will be returned
/// @param location  An optional string to indicate where this value is used
template<typename T>
T getNoBoundary(const Options& option, const std::string& location = "") {
#if CHECKLEVEL >= 1
  // Mark option as final inside the domain
  const_cast<Options&>(option).attributes["final-domain"] = location;
#endif
  return getNonFinal<T>(option);
}

#if CHECKLEVEL >= 1
/// A wrapper around get<>() which captures debugging information
///
/// Usage:
///   auto var = GET_NOBOUNDARY(Field3D, option["value"]);
#define GET_NOBOUNDARY(Type, option) \
  getNoBoundary<Type>(option, __FILE__ ":" TOSTRING(__LINE__))
#else
#define GET_NOBOUNDARY(Type, option) \
  getNoBoundary<Type>(option)
#endif

/// Set values in an option. This could be optimised, but
/// currently the is_value private variable would need to be modified.
///
/// If the value has been used then raise an exception (if CHECK >= 1)
/// This is to prevent values being modified after use.
///
/// @tparam T The type of the value to set. Usually this is inferred
template<typename T>
Options& set(Options& option, T value) {
  // Check that the value has not already been used
#if CHECKLEVEL >= 1
  if (option.hasAttribute("final")) {
    throw BoutException("Setting value of {} but it has already been used in {}.",
                        option.name(), option.attributes["final"].as<std::string>());
  }
  if (option.hasAttribute("final-domain")) {
    throw BoutException("Setting value of {} but it has already been used in {}.",
                        option.name(),
                        option.attributes["final-domain"].as<std::string>());
  }

#endif
  option.force(std::move(value));
  return option;
}

/// Set values in an option. This could be optimised, but
/// currently the is_value private variable would need to be modified.
///
/// This version only checks that the boundary cells have not
/// already been used by a call to get, not a call to getNoBoundary
/// or getNonFinal.
///
/// @tparam T The type of the value to set. Usually this is inferred
template<typename T>
Options& setBoundary(Options& option, T value) {
  // Check that the value has not already been used
#if CHECKLEVEL >= 1
  if (option.hasAttribute("final")) {
    throw BoutException("Setting boundary of {} but it has already been used in {}.",
                        option.name(), option.attributes["final"].as<std::string>());
  }
#endif
  option.force(std::move(value));
  return option;
}

/// Add value to a given option. If not already set, treats
/// as zero and sets the option to the value.
///
/// @tparam T The type of the value to add. The existing value
///           will be casted to this type
///
/// @param option  The value to modify (or set if not already set)
/// @param value   The quantity to add.
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
///
/// @param option  The value to modify (or set if not already set)
/// @param value   The quantity to add.
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

template<typename T>
void set_with_attrs(Options& option, T value, std::initializer_list<std::pair<std::string, Options::AttributeType>> attrs) {
  option = value;
  option.setAttributes(attrs);
}

#endif // HERMES_COMPONENT_H
