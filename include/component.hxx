#pragma once

#ifndef HERMES_COMPONENT_H
#define HERMES_COMPONENT_H

#include <cmath>
#include <initializer_list>
#include <map>
#include <memory>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

#include <bout/assert.hxx>
#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/field2d.hxx>
#include <bout/field3d.hxx>
#include <bout/generic_factory.hxx>
#include <bout/options.hxx>
#include <bout/unused.hxx>
#include <fmt/args.h>

#include "guarded_options.hxx"
#include "permissions.hxx"
#include "hermes_utils.hxx"

class Solver; // Time integrator

/// Simple struct to store information on the different types of
/// species present in a simulation
struct SpeciesInformation {
  SpeciesInformation(const std::vector<std::string>& electrons,
                     const std::vector<std::string>& neutrals,
                     const std::vector<std::string>& positive_ions,
                     const std::vector<std::string> & negative_ions)
    : electrons(electrons), neutrals(neutrals), positive_ions(positive_ions), negative_ions(negative_ions), ions(positive_ions) {
    finish_construction();
  }

  SpeciesInformation(const std::initializer_list<std::string> species) {
    for (auto& sp : species) {
      // FIXME: identifySpecies only identifies positive ions
      // FIXME: identifySpecies has no concept of ebeam
      const SpeciesType type = identifySpeciesType(sp);
      if (type == SpeciesType::electron) {
        electrons.push_back(sp);
      } else if (type == SpeciesType::ion) {
        positive_ions.push_back(sp);
      } else if (type == SpeciesType::neutral) {
        neutrals.push_back(sp);
      } else {
        throw BoutException("Species {} has unrecognised type {}", sp, toString(type));
      }
      finish_construction();
    }
  }

  std::vector<std::string> electrons, neutrals, positive_ions, negative_ions, ions, charged, non_electrons, all_species;

  private:
    void finish_construction() {
      ions = positive_ions;
      ions.insert(ions.end(), negative_ions.begin(), negative_ions.end());
      charged = ions;
      charged.insert(charged.end(), electrons.begin(), electrons.end());
      non_electrons = ions;
      non_electrons.insert(non_electrons.end(), neutrals.begin(), neutrals.end());
      all_species = charged;
      all_species.insert(all_species.end(), neutrals.begin(), neutrals.end());
    }
};

/// Interface for a component of a simulation model
///
/// The constructor of derived types should have signature
///   (std::string name, Options &options, Solver *solver)
///
struct Component {
  /// Initialise the `state_variable_acceess` permissions. Note that
  /// `{all_species}` in any variable names will be replaced with the
  /// names of all species being simulated (by calling
  /// `declareAllSpecies()`, which is done after all components are
  /// created by a ComponentSchedular).
  Component(Permissions&& access_permissions)
      : state_variable_access(access_permissions) {}

  virtual ~Component() {}

  /// Modify the given simulation state. This method will wrap the
  /// state in a GuardedOptions object and pass that to the private
  /// implementation of transform provided by each component.
  void transform(Options &state);
  
  /// Use the final simulation state to update internal state
  /// (e.g. time derivatives)
  virtual void finally(const Options &UNUSED(state)) { }

  /// Add extra fields for output, or set attributes e.g docstrings
  virtual void outputVars(Options &UNUSED(state)) { }

  /// Add extra fields to restart files
  virtual void restartVars(Options &UNUSED(state)) { }

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

  /// Perform all the substitutions registered with this component. In addition,
  /// the following substitutions will be done using the information on the
  /// species names in the simulation:
  ///   - electrons (any electron species)
  ///   - electrons2 (same as above, used for cross-product)
  ///   - neutrals (species with no charge)
  ///   - neutrals2 (same as above, used for cross-product)
  ///   - positive_ions (ions with a positive charge)
  ///   - positive_ions2 (same as above, used for cross-product)
  ///   - negative_ions (ions with a negative charge)
  ///   - negative_ions2 (same as above, used for cross-product)
  ///   - ions (all ions, regardless of sign of charge)
  ///   - ions2 (same as above, used for cross-product)
  ///   - charged (ions and electrons)
  ///   - charged2 (same as above, used for cross-product)
  ///   - non_electrons (ions and neutrals)
  ///   - non_electrons2 (same as above, used for cross-product)
  ///   - all_species (ions, neutrals, and electrons)
  ///   - all_species2 (same as above, used for cross-product)
  ///
  /// Note that substitutions must be provied for *all* labels present
  /// in the permissions data, or else an exception will be thrown.
  void performAllSubstitutions(const SpeciesInformation& info);

protected:
  /// Set the level of access needed by this component for a particular variable.
  void setPermissions(const std::string& variable,
                      const Permissions::AccessRights& rights) {
    state_variable_access.setAccess(variable, rights);
  }
  void setPermissions(const Permissions::VarRights& info) {
    setPermissions(info.name, info.rights);
  }

  /// Replace a substitutions to be performed on the placeholder in
  /// the name of variables stored in the access control information
  /// for this component. The substitution won't actually happen until
  /// `performAllSubstitutions` is called. If you register the same
  /// label more then once then an exception will be
  /// thrown. Therefore, be careful not to use any of the built-in
  /// labels that will be set by `performAllSubstitutions`.
  void substitutePermissions(const std::string& label,
                             const std::vector<std::string>& substitution) {
    registerSubstitution(label, substitution);
  }

private:
  /// Information on which state variables the transform method will read and write.
  Permissions state_variable_access;
  std::map<std::string, std::vector<std::string>> permission_substitutions;

  /// Modify the given simulation state. All components must
  /// implement this function. It will only allow the reading
  /// from/writing to state variables with the appropriate permissiosn
  /// in `state_variable_access`.
  virtual void transform_impl(GuardedOptions &state) = 0;

  void registerSubstitution(const std::string& label,
                            const std::vector<std::string>& substitution) {
    if (permission_substitutions.count(label) > 0) {
      throw BoutException("Substitution with label {} already present", label);
    }
    permission_substitutions[label] = substitution;
  }
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
template<typename T>
T getNonFinal(const GuardedOptions & option) {
  return getNonFinal<T>(option.get());
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
template<typename T>
T get(const GuardedOptions & option, const std::string& location = "") {
  return get<T>(option.get(), location);
}

/// Check if an option can be fetched
/// Sets the final flag so setting the value
/// afterwards will lead to an error
bool isSetFinal(const Options& option, const std::string& location = "");
bool isSetFinal(const GuardedOptions & option, const std::string& location = "");

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
bool isSetFinalNoBoundary(const GuardedOptions & option, const std::string& location = "");

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
template<typename T>
T getNoBoundary(const GuardedOptions & option, const std::string& location = "") {
  return getNoBoundary<T>(option.get(Regions::Interior), location);
}
template<typename T>
T getNoBoundary(const GuardedOptions && option, const std::string& location = "") {
  return getNoBoundary<T>(option.get(Regions::Interior), location);
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

/// Check whether value is valid, returning true
/// if invalid i.e contains non-finite values
template <typename T>
bool hermesDataInvalid([[maybe_unused]] const T& value) {
  return false; // Default
}

/// Check Field3D values.
/// Doesn't check boundary cells
template<>
inline bool hermesDataInvalid(const Field3D& value) {
  for (auto& i : value.getRegion("RGN_NOBNDRY")) {
    if (!std::isfinite(value[i])) {
      return true;
    }
  }
  return false;
}

/// Check Field2D values.
/// Doesn't check boundary cells
template<>
inline bool hermesDataInvalid(const Field2D& value) {
  for (auto& i : value.getRegion("RGN_NOBNDRY")) {
    if (!std::isfinite(value[i])) {
      return true;
    }
  }
  return false;
}

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

  if (hermesDataInvalid(value)) {
    throw BoutException("Setting invalid value for '{}'", option.str());
  }
#endif

  option.force(std::move(value));
  return option;
}
template<typename T>
GuardedOptions & set(GuardedOptions & option, T value) {
  set(option.getWritable(), value);
  return option;
}
template<typename T>
GuardedOptions && set(GuardedOptions && option, T value) {
  set(option.getWritable(), value);
  return std::move(option);
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
template<typename T>
GuardedOptions & setBoundary(GuardedOptions & option, T value) {
  setBoundary(option.getWritable(Regions::Boundaries), value);
  return option;
}
template<typename T>
GuardedOptions && setBoundary(GuardedOptions && option, T value) {
  setBoundary(option.getWritable(Regions::Boundaries), value);
  return std::move(option);
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
template<typename T>
GuardedOptions & add(GuardedOptions & option, T value) {
  add(option.getWritable(), value);
  return option;
}
template<typename T>
GuardedOptions && add(GuardedOptions && option, T value) {
  add(option.getWritable(), value);
  return std::move(option);
}

/// Add value to a given option. If not already set, treats
/// as zero and sets the option to the value.
///
/// @param option  The value to modify (or set if not already set)
/// @param value   The quantity to add.
template<typename T>
Options& subtract(Options& option, T value) {
  if (!option.isSet()) {
    return set(option, -value);
  } else {
    try {
      return set(option, bout::utils::variantStaticCastOrThrow<Options::ValueType, T>(option.value) - value);
    } catch (const std::bad_cast &e) {
      // Convert to a more useful error message
      throw BoutException("Could not convert {:s} to type {:s}",
                          option.str(), typeid(T).name());
    }
  }
}
template<typename T>
GuardedOptions & subtract(GuardedOptions & option, T value) {
  subtract(option.getWritable(), value);
  return option;
}
template<typename T>
GuardedOptions && subtract(GuardedOptions && option, T value) {
  subtract(option.getWritable(), value);
  return std::move(option);
}

template<typename T>
void set_with_attrs(Options& option, T value, std::initializer_list<std::pair<std::string, Options::AttributeType>> attrs) {
  option.force(value);
  option.setAttributes(attrs);
}
template<typename T>
void set_with_attrs(GuardedOptions & option, T value, std::initializer_list<std::pair<std::string, Options::AttributeType>> attrs) {
  set_with_attrs(option.getWritable(), value, attrs);
}
template<typename T>
void set_with_attrs(GuardedOptions && option, T value, std::initializer_list<std::pair<std::string, Options::AttributeType>> attrs) {
  set_with_attrs(option.getWritable(), value, attrs);
}

#if CHECKLEVEL >= 1
template<>
inline void set_with_attrs(Options& option, Field3D value, std::initializer_list<std::pair<std::string, Options::AttributeType>> attrs) {
  if (!value.isAllocated()) {
    throw BoutException("set_with_attrs: Field3D assigned to {:s} is not allocated", option.str());
  }
  option.force(value);
  option.setAttributes(attrs);
}
template<>
inline void set_with_attrs(GuardedOptions & option, Field3D value, std::initializer_list<std::pair<std::string, Options::AttributeType>> attrs) {
  set_with_attrs(option.getWritable(), std::move(value), attrs);
}
template<>
inline void set_with_attrs(GuardedOptions && option, Field3D value, std::initializer_list<std::pair<std::string, Options::AttributeType>> attrs) {
  set_with_attrs(option.getWritable(), std::move(value), attrs);
}
#endif

#endif // HERMES_COMPONENT_H
