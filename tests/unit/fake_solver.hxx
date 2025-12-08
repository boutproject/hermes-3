#pragma once

#include <bout/solver.hxx>

class FakeSolver : public Solver {
public:
  /// Set the state variables
  void setState(Options state) {
    for (auto& v : f3d) {
      *(v.var) = state[v.name].as<Field3D>();
    }
  }

  Options getState() const {
    Options state;
    for (const auto& v : f3d) {
      state[v.name] = *(v.var);
    }
    return state;
  }

  /// Get the time derivative fields
  Options getTimeDerivs() const {
    Options ddts;
    for (const auto& v : f3d) {
      ddts[v.name] = *(v.F_var);
    }
    return ddts;
  }

  int run() override { return 1; }
};
