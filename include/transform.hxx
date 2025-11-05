#pragma once
#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "component.hxx"

/// Apply changes to the state
///
struct Transform : public Component {
  Transform(std::string name, Options& options, Solver*);

private:
  std::map<std::string, std::string> transforms;

  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<Transform> registercomponenttransform("transform");
}

#endif // TRANSFORM_H
