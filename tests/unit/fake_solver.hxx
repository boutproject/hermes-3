#pragma once

#include <bout/solver.hxx>
#include <vector>

class FakeSolver : public Solver {
public:
  bool evolvesBoundary(const std::string& name) const {
    for (const auto& v : f3d) {
      if (v.name == name) {
        return v.evolve_bndry;
      }
    }
    return false;
  }

  void applyValueBoundaries(BoutReal time) {
    for (const auto& v : f3d) {
      v.var->applyBoundary(time);
    }
  }
  void applyDerivativeBoundaries() {
    for (const auto& v : f3d) {
      if (v.evolve_bndry) {
        v.var->applyTDerivBoundary();
      }
    }
  }

  std::vector<BoutReal> packState() {
    size_t count = 0;
    for (const auto& v : f3d) {
      count += size(v.var->getRegion("RGN_NOBNDRY"));
      if (v.evolve_bndry) {
        count += size(v.var->getRegion("RGN_BNDRY"));
      }
    }
    std::vector<BoutReal> packed(count);
    save_vars(packed.data());
    return packed;
  }
  void unpackState(std::vector<BoutReal>& packed) { load_vars(packed.data()); }

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
    for (const auto& v : f2d) {
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
    for (const auto& v : f2d) {
      ddts[v.name] = *(v.F_var);
    }
    return ddts;
  }

  int run() override { return 1; }
};
