#include <memory>
#include <string>

#include <bout/assert.hxx>
#include <bout/options.hxx>
#include <bout/output.hxx>

#include "../include/component.hxx"
#include "../include/guarded_options.hxx"
#include "../include/permissions.hxx"

std::unique_ptr<Component> Component::create(const std::string& type,
                                             const std::string& name, Options& alloptions,
                                             Solver* solver) {

  return ComponentFactory::getInstance().create(type, name, alloptions, solver);
}

void Component::transform(Options& state) {
  GuardedOptions guarded(&state, &state_variable_access);
  transform_impl(guarded);
#if CHECKLEVEL >= 999
  for (auto& [varname, region] : guarded.unreadItems()) {
    output_warn.write("Did not read from state variable {} in region(s) {}\n", varname,
                      Permissions::regionNames(region));
  }
  for (auto& [varname, region] : guarded.unwrittenItems()) {
    output_warn.write("Did not write to state variable {} in region(s) {}\n", varname,
                      Permissions::regionNames(region));
  }
#endif
}

void Component::declareAllSpecies(const SpeciesInformation& info) {
  state_variable_access.substitute("electrons", info.electrons);
  state_variable_access.substitute("electrons2", info.electrons);
  state_variable_access.substitute("neutrals", info.neutrals);
  state_variable_access.substitute("neutrals2", info.neutrals);
  state_variable_access.substitute("positive_ions", info.positive_ions);
  state_variable_access.substitute("positive_ions2", info.positive_ions);
  state_variable_access.substitute("negative_ions", info.negative_ions);
  state_variable_access.substitute("negative_ions2", info.negative_ions);
  state_variable_access.substitute("ions", info.ions);
  state_variable_access.substitute("ions2", info.ions);
  state_variable_access.substitute("charged", info.charged);
  state_variable_access.substitute("charged2", info.charged);
  state_variable_access.substitute("non_electrons", info.non_electrons);
  state_variable_access.substitute("non_electrons2", info.non_electrons);
  state_variable_access.substitute("all_species", info.all_species);
  state_variable_access.substitute("all_species2", info.all_species);
  state_variable_access.checkNoRemainingSubstitutions();
}

constexpr decltype(ComponentFactory::type_name) ComponentFactory::type_name;
constexpr decltype(ComponentFactory::section_name) ComponentFactory::section_name;
constexpr decltype(ComponentFactory::option_name) ComponentFactory::option_name;
constexpr decltype(ComponentFactory::default_type) ComponentFactory::default_type;

bool isSetFinal(const Options& option, [[maybe_unused]] const std::string& location) {
#if CHECKLEVEL >= 1
  // Mark option as final, both inside the domain and the boundary
  const_cast<Options&>(option).attributes["final"] = location;
  const_cast<Options&>(option).attributes["final-domain"] = location;
#endif
  return option.isSet();
}

bool isSetFinal(const GuardedOptions& option,
                [[maybe_unused]] const std::string& location) {
  const bool set = option.isSet();
#if CHECKLEVEL >= 1
  const PermissionTypes perm = option.getHighestPermission();
  if (perm >= PermissionTypes::Read or (perm == PermissionTypes::ReadIfSet and set)) {
    const Options& opt = option.get();
    const_cast<Options&>(opt).attributes["final"] = location;
    const_cast<Options&>(opt).attributes["final-domain"] = location;
  }
#endif
  return set;
}

bool isSetFinalNoBoundary(const Options& option,
                          [[maybe_unused]] const std::string& location) {
#if CHECKLEVEL >= 1
  // Mark option as final inside the domain, but not in the boundary
  const_cast<Options&>(option).attributes["final-domain"] = location;
#endif
  return option.isSet();
}

bool isSetFinalNoBoundary(const GuardedOptions& option,
                          [[maybe_unused]] const std::string& location) {
  const bool set = option.isSet();
#if CHECKLEVEL >= 1
  const PermissionTypes perm = option.getHighestPermission(Regions::Interior);
  if (perm >= PermissionTypes::Read or (perm == PermissionTypes::ReadIfSet and set)) {
    // Mark option as final inside the domain, but not in the boundary
    const_cast<Options&>(option.get(Regions::Interior)).attributes["final-domain"] =
        location;
  }
#endif
  return set;
}

/// Formatter for Components
///
/// By default it will use the format `OBJECT_NAME
/// (COMPONENT_TYPE_NAME)`, if the two names are different. If they
/// are the same then it will just show the object name. This
/// behaviour can be overriden using the format specifiers below:
///
/// - ``~n``: Don't show object name
/// - ``~t``: Don't show type name
/// - ``T``: Always show type name
constexpr auto fmt::formatter<Component>::parse(format_parse_context& ctx)
    -> format_parse_context::iterator {
  const auto* it = ctx.begin();
  const auto* end = ctx.end();

  if (it == end) {
    return underlying.parse(ctx);
  }

  while (it != end and *it != ':' and *it != '}') {
    // Other cases handled explicitly below
    // NOLINTNEXTLINE(bugprone-switch-missing-default-case)
    switch (*it) {
    case '~':
      ++it;
      if (*it == 'n') {
        hide_name = true;
      } else if (*it == 't') {
        hide_type = true;
      } else {
        throw fmt::format_error(fmt::format("Unrecognised format specifier '~{}'", *it));
      }
      ++it;
      break;
    case 'T':
      show_type = true;
      ++it;
      break;
    }
    // FIXME: This is copied from the parser used for Fields, but
    // won't it result in an infinite loop if an unrecognised
    // character is encountered?
  }

  if (hide_type and show_type) {
    throw fmt::format_error("Format specifiers 'T' and '~t' are mutually-exclusive");
  }

  if (it != end and *it != '}') {
    if (*it != ':') {
      throw fmt::format_error("invalid format specifier");
    }
    ++it;
  }

  ctx.advance_to(it);
  return underlying.parse(ctx);
}

auto fmt::formatter<Component>::format(const Component& component,
                                       format_context& ctx) const
    -> format_context::iterator {
  const std::string tn = component.typeName();
  const std::string on = component.objectName();

  std::string result;

  if (not hide_name) {
    result += on;
  }
  if (show_type or (tn != on and not hide_type)) {
    if (not hide_name) {
      result += " (";
    }
    result += tn;
    if (not hide_name) {
      result += ")";
    }
  }
  return underlying.format(result, ctx);
}
