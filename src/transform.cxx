
#include "../include/transform.hxx"

#include <bout/utils.hxx> // for trim, strsplit

Transform::Transform(std::string name, Options& alloptions, Solver* UNUSED(solver)) {

  Options& options = alloptions[name];

  const auto trim_chars = " \t\r()";

  const auto str = trim(
      options["transforms"].doc("Comma-separated list e.g. a = b, c = d"), trim_chars);

  for (const auto& assign_str : strsplit(str, ',')) {
    auto assign_lr = strsplit(assign_str, '=');
    if (assign_lr.size() != 2) {
      throw BoutException("Expected one assignment ('=') in '{}'", assign_str);
    }

    const auto left = trim(assign_lr.front(), trim_chars);
    const auto right = trim(assign_lr.back(), trim_chars);

    transforms[left] = right;
  }
}

void Transform::transform(Options& state) {
  for (const auto& lr : transforms) {
    state[lr.first] = state[lr.second];
  }
}
