#pragma once
#ifndef NOFLOW_BOUNDARY_H
#define NOFLOW_BOUNDARY_H

#include "component.hxx"

struct NoFlowBoundary : public NamedComponent<NoFlowBoundary> {
  NoFlowBoundary(std::string name, Options& alloptions, Solver*)
      : NamedComponent(name, {writeBoundaryIfSet("species:{name}:{variables}")}) {

    Options& options = alloptions[name];
    noflow_lower_y = options["noflow_lower_y"]
                         .doc("No-flow boundary on lower y?")
                         .withDefault<bool>(true);
    noflow_upper_y = options["noflow_upper_y"]
                         .doc("No-flow boundary on upper y?")
                         .withDefault<bool>(true);
    substitutePermissions("name", {name});
    substitutePermissions("variables",
                          {"density", "temperature", "pressure", "velocity", "momentum"});
  }

  static constexpr auto type = "noflow_boundary";

private:
  std::string name;    ///<
  bool noflow_lower_y; ///< No-flow boundary on lower y?
  bool noflow_upper_y; ///< No-flow boundary on upper y?

  /// Inputs
  ///  - species
  ///    - <name>
  ///      - density      [Optional]
  ///      - temperature  [Optional]
  ///      - pressure     [Optional]
  ///      - velocity     [Optional]
  ///      - momentum     [Optional]
  void transform_impl(GuardedOptions& state) override;
};

namespace {
RegisterComponent<NoFlowBoundary> registercomponentnoflowboundary;
}

#endif // NOFLOW_BOUNDARY_H
