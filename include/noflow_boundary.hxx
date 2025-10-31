#pragma once
#ifndef NOFLOW_BOUNDARY_H
#define NOFLOW_BOUNDARY_H

#include "component.hxx"

struct NoFlowBoundary : public Component {
  NoFlowBoundary(std::string name, Options& alloptions, Solver*) : name(name) {
    AUTO_TRACE();

    Options& options = alloptions[name];
    noflow_lower_y = options["noflow_lower_y"]
                         .doc("No-flow boundary on lower y?")
                         .withDefault<bool>(true);
    noflow_upper_y = options["noflow_upper_y"]
                         .doc("No-flow boundary on upper y?")
                         .withDefault<bool>(true);
  }

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
  void transform(GuardedOptions& state) override;
};

namespace {
RegisterComponent<NoFlowBoundary> registercomponentnoflowboundary("noflow_boundary");
}

#endif // NOFLOW_BOUNDARY_H
