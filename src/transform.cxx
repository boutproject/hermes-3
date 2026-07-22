
#include "../include/transform.hxx"

#include <bout/utils.hxx> // for trim, strsplit

Transform::Transform(std::string name, Options& alloptions, Solver* UNUSED(solver))
    : NamedComponent(name, {readOnly("{inputs}"), writeFinal("{outputs}")}) {

  Options& options = alloptions[name];

  const auto trim_chars = " \t\r()";

  const auto str = trim(
      options["transforms"].doc("Comma-separated list e.g. a = b, c = d"), trim_chars);

  std::vector<std::string> inputs, outputs;

  for (const auto& assign_str : strsplit(str, ',')) {
    auto assign_lr = strsplit(assign_str, '=');
    if (assign_lr.size() != 2) {
      throw BoutException("Expected one assignment ('=') in '{}'", assign_str);
    }

    const auto left = trim(assign_lr.front(), trim_chars);
    const auto right = trim(assign_lr.back(), trim_chars);

    transforms[left] = right;
    inputs.push_back(right);
    outputs.push_back(left);
  }

  substitutePermissions("inputs", inputs);
  substitutePermissions("outputs", outputs);
}

void Transform::transform_impl(GuardedOptions& state) {
  for (const auto& lr : transforms) {
    // FIXME: The assignment of one Options object to another means
    // that the name of the first one is copied. This is a problem
    // because the name now checked against permissions can differ
    // from the name used to access the data. Ideally Options would
    // rewrite the name on assignment, but possibly there's a reason
    // it does't?
    //
    // Choices:
    // - Rewrite assignment for Options so it rewrites names if one is already present.
    // - Add a setter method for full_name and use this to correct
    //   - The correction could be done in an assignment operator for GuardedOptions or else here
    // - Copy only the underlying object, rather than the option itself. This will bake in an assumption about the type though.
    state[lr.first].getWritable() = state[lr.second].get().copy();
  }
}
