#pragma once
#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "component.hxx"

/// Apply changes to the state
///
struct Transform : public Component {
  Transform(std::string name, Options& options, Solver*);

  void transform(Options& state) override;

private:
  std::map<std::string, std::string> transforms;
};

namespace {
RegisterComponent<Transform> registercomponenttransform("transform");
}

#endif // TRANSFORM_H
