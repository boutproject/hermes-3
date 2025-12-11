#pragma once
#ifndef HERMES_COMPONENT_H
#define HERMES_COMPONENT_H

#include <cmath>
#include <initializer_list>
#include <map>
#include <memory>
#include <string>
#include <type_traits>
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
                     const std::vector<std::string>& negative_ions)
      : electrons(std::move(electrons)), neutrals(std::move(neutrals)),
        positive_ions(std::move(positive_ions)), negative_ions(std::move(negative_ions)),
        ions(std::move(positive_ions)) {
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

/// Lightweight type to store information used to build a component.
struct ComponentInformation {
  std::string name;
  std::string type;

  ComponentInformation(const std::string& name_, const std::string& type_) : name(name_), type(type_) {}

  ComponentInformation(std::string&& name_, std::string&& type_) : name(std::move(name_)), type(std::move(type_)) {}

  bool operator<(const ComponentInformation& other) const {
    return std::pair(name, type) < std::pair(other.name, other.type);
  }

  bool operator==(const ComponentInformation& other) const {
    return std::pair(name, type) == std::pair(other.name, other.type);
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

  /// Return a list of names/types of other components needed by this
  /// component. All configurations for these components will take the
  /// default value, unless set in the input file.
  virtual std::vector<ComponentInformation> additionalComponents() { return {}; }

  /// Modify the given simulation state. This method will wrap the
  /// state in a GuardedOptions object and pass that to the private
  /// implementation of transform provided by each component.
  void transform(Options& state);

  /// Use the final simulation state to update internal state
  /// (e.g. time derivatives)
  virtual void finally(const Options& UNUSED(state)) {}

  /// Add extra fields for output, or set attributes e.g docstrings
  virtual void outputVars(Options& UNUSED(state)) {}

  /// Add extra fields to restart files
  virtual void restartVars(Options& UNUSED(state)) {}

  /// Preconditioning
  virtual void precon(const Options& UNUSED(state), BoutReal UNUSED(gamma)) {}

  /// Create a Component
  ///
  /// @param type     The name of the component type to create (e.g. "evolve_density")
  /// @param name     The species/name for this instance.
  /// @param options  Component settings: options[name] are specific to this component
  /// @param solver   Time-integration solver
  static std::unique_ptr<Component> create(const std::string& type,
                                           const std::string& name, Options& options,
                                           Solver* solver);

  /// Tell the component the name of all species in the simulation, by type. It
  /// will use this information to substitute the following placeholders in
  /// `state_variable_access`:
  ///   - electrons (any electron species)
  ///   - electrons2 (same as above, used for Cartesian product)
  ///   - neutrals (species with no charge)
  ///   - neutrals2 (same as above, used for Cartesian product)
  ///   - positive_ions (ions with a positive charge)
  ///   - positive_ions2 (same as above, used for Cartesian product)
  ///   - negative_ions (ions with a negative charge)
  ///   - negative_ions2 (same as above, used for Cartesian product)
  ///   - ions (all ions, regardless of sign of charge)
  ///   - ions2 (same as above, used for Cartesian product)
  ///   - charged (ions and electrons)
  ///   - charged2 (same as above, used for Cartesian product)
  ///   - non_electrons (ions and neutrals)
  ///   - non_electrons2 (same as above, used for Cartesian product)
  ///   - all_species (ions, neutrals, and electrons)
  ///   - all_species2 (same as above, used for Cartesian product)
  ///
  /// At the end of this function there is a call to
  /// Permissions::checkNoRemainingSubstitutions. All substitutions
  /// must be completed or else an exception will be thrown.
  void declareAllSpecies(const SpeciesInformation & info);

  const Permissions& getPermissions() const { return state_variable_access; }

protected:
  /// Set the level of access needed by this component for a particular variable.
  void setPermissions(const std::string& variable,
                      const Permissions::AccessRights& rights) {
    state_variable_access.setAccess(variable, rights);
  }
  void setPermissions(const Permissions::VarRights& info) {
    setPermissions(info.name, info.rights);
  }

  /// Replace a placeholder in the name of variables stored in the access control
  /// information for this component.
  void substitutePermissions(const std::string& label,
                             const std::vector<std::string>& substitution) {
    state_variable_access.substitute(label, substitution);
  }

private:
  /// Information on which state variables the transform method will read and write.
  Permissions state_variable_access;

  /// The implementation of the transform method. Modify the given
  /// simulation state. All components must implement this
  /// function. It will only allow the reading from/writing to state
  /// variables with the appropriate permissiosn in
  /// `state_variable_access`.
  virtual void transform_impl(GuardedOptions &state) = 0;
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


namespace hermes {
/// Enable a function if and only if `T` is a (subclass of) `GuardedOptions`
template <class T>
using EnableIfGuardedOption = std::enable_if_t<std::is_base_of_v<GuardedOptions, T>>;
}

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
T get(const Options& option, [[maybe_unused]] const std::string& location = "") {
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
T getNoBoundary(const Options& option, [[maybe_unused]] const std::string& location = "") {
#if CHECKLEVEL >= 1
  // Mark option as final inside the domain
  const_cast<Options&>(option).attributes["final-domain"] = location;
#endif
  return getNonFinal<T>(option);
}

template<typename T, class GO, typename = hermes::EnableIfGuardedOption<GO>>
T getNoBoundary(GO&& option, const std::string& location = "") {
  return getNoBoundary<T>(std::forward<GO>(option).get(Regions::Interior), location);
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

template<typename T, class GO, typename = hermes::EnableIfGuardedOption<GO>>
decltype(auto) set(GO&& option, T value) {
  set(std::forward<GO>(option).getWritable(), value);
  return std::forward<GO>(option);
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

template<typename T, class GO, typename = hermes::EnableIfGuardedOption<GO>>
decltype(auto) setBoundary(GO&& option, T value) {
  setBoundary(std::forward<GO>(option).getWritable(Regions::Boundaries), value);
  return std::forward<GO>(option);
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
  }
  try {
    return set(option, value
                           + bout::utils::variantStaticCastOrThrow<Options::ValueType, T>(
                               option.value));
  } catch (const std::bad_cast& e) {
    // Convert to a more useful error message
    throw BoutException("Could not convert {:s} to type {:s}", option.str(),
                        typeid(T).name());
  }
}

template<typename T, class GO, typename = hermes::EnableIfGuardedOption<GO>>
decltype(auto) add(GO&& option, T value) {
  add(std::forward<GO>(option).getWritable(), value);
  return std::forward<GO>(option);
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
  }
  try {
    return set(option,
               bout::utils::variantStaticCastOrThrow<Options::ValueType, T>(option.value)
                   - value);
  } catch (const std::bad_cast& e) {
    // Convert to a more useful error message
    throw BoutException("Could not convert {:s} to type {:s}", option.str(),
                        typeid(T).name());
  }
}

template<typename T, class GO, typename = hermes::EnableIfGuardedOption<GO>>
decltype(auto) subtract(GO&& option, T value) {
  subtract(std::forward<GO>(option).getWritable(), value);
  return std::forward<GO>(option);
}

template<typename T>
void set_with_attrs(Options& option, T value, std::initializer_list<std::pair<std::string, Options::AttributeType>> attrs) {
  option.force(value);
  option.setAttributes(attrs);
}

template<typename T, class GO, typename = hermes::EnableIfGuardedOption<GO>>
void set_with_attrs(GO&& option, T value, std::initializer_list<std::pair<std::string, Options::AttributeType>> attrs) {
  set_with_attrs(std::forward<GO>(option).getWritable(), value, attrs);
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

template<class GO, typename = hermes::EnableIfGuardedOption<GO>>
inline void set_with_attrs(GO&& option, Field3D value, std::initializer_list<std::pair<std::string, Options::AttributeType>> attrs) {
  set_with_attrs(std::forward<GO>(option).getWritable(), std::move(value), attrs);
}
#endif

#endif // HERMES_COMPONENT_H
