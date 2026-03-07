#pragma once
#ifndef NOFLOW_BOUNDARY_H
#define NOFLOW_BOUNDARY_H

#include "component.hxx"

struct NoFlowBoundary : public Component {
  NoFlowBoundary(std::string name, Options& alloptions, Solver*)
      : Component({writeBoundaryIfSet("species:{name}:{variables}")}), name(name) {

    Options& options = alloptions[name];
    noflow_lower_y = options["noflow_lower_y"]
                         .doc("No-flow boundary on lower y?")
                         .withDefault<bool>(true);
    noflow_upper_y = options["noflow_upper_y"]
                         .doc("No-flow boundary on upper y?")
                         .withDefault<bool>(true);
    noflow_inner_x = options["noflow_lower_y"]
                         .doc("No-flow boundary on inner x?")
                         .withDefault<bool>(false);
    noflow_outer_x = options["noflow_upper_y"]
                         .doc("No-flow boundary on outer x?")
                         .withDefault<bool>(false);

    substitutePermissions("name", {name});
    substitutePermissions("variables",
                          {"density", "temperature", "pressure", "velocity", "momentum"});
  }

private:
  std::string name;    ///<
  bool noflow_lower_y; ///< No-flow boundary on lower y?
  bool noflow_upper_y; ///< No-flow boundary on upper y?
  bool noflow_inner_x; ///< No-flow boundary on inner x?
  bool noflow_outer_x; ///< No-flow boundary on outer x?

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
RegisterComponent<NoFlowBoundary> registercomponentnoflowboundary("noflow_boundary");
}

#endif // NOFLOW_BOUNDARY_H
