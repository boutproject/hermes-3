
#include "../include/transform.hxx"

#include <bout/utils.hxx> // for trim, strsplit

Transform::Transform(std::string name, Options& alloptions, Solver* UNUSED(solver))
    : Component({readOnly("{inputs}"), writeFinal("{outputs}")}) {

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
    inputs.push_back(left);
    outputs.push_back(right);
  }

  state_variable_access.substitute("inputs", inputs);
  state_variable_access.substitute("outputs", outputs);
}

void Transform::transform_impl(GuardedOptions& state) {
  for (const auto& lr : transforms) {
    state[lr.first].getWritable() = state[lr.second].get().copy();
  }
}
