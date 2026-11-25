#pragma once

#include "component.hxx"

#include <bout/boutcomm.hxx>
#include <bout/constants.hxx>
#include <bout/mesh.hxx>
#include <bout/options.hxx>
#include <bout/output_bout_types.hxx>

#include <map>
#include <string>
#include <vector>

/// Small, MPI-reduced ledgers of the terms actually used in the equations.
/// Integrals exclude guards; conversion includes the simulated toroidal extent.
/// No flux operator is evaluated here. Collectives run only at output times.
struct DomainBalance {
  std::map<std::pair<std::string, std::string>, BoutReal> local;

  void clear() { local.clear(); }

  static BoutReal integral(const Field3D& field) {
    ASSERT1(field.isAllocated());
    const auto* c = field.getCoordinates();
    BoutReal result = 0.0;
    for (const auto& i : field.getRegion("RGN_NOBNDRY")) {
      result += field[i] * c->J()[i] * c->dx()[i] * c->dy()[i] * c->dz()[i];
    }
    return result;
  }

  void set(const std::string& quantity, const std::string& term, const Field3D& field) {
    local[{quantity, term}] = integral(field);
  }

  void add(const std::string& quantity, const std::string& term, BoutReal value) {
    local[{quantity, term}] += value;
  }

  static std::string target(Mesh* mesh, bool lower) {
    if (lower) {
      return mesh->getGlobalYIndexNoBoundaries(mesh->ystart) == 0 ? "lower_y"
                                                                  : "lower_y_inner";
    }
    const int last = mesh->GlobalNy - (mesh->LocalNy - mesh->yend);
    return mesh->getGlobalYIndex(mesh->yend) == last ? "upper_y" : "upper_y_inner";
  }

  void initTargets(const std::string& quantity, const std::string& prefix) {
    for (const auto* side : {"lower_y", "upper_y", "lower_y_inner", "upper_y_inner"}) {
      local.try_emplace({quantity, prefix + "_" + side}, 0.0);
      local.try_emplace({quantity, "present_" + std::string{side}}, 0.0);
    }
  }

  /// flow is the production operator's area-integrated lower-face flow.
  void addTargets(const std::string& quantity, const std::string& prefix,
                  const Field3D& flow) {
    initTargets(quantity, prefix);
    if (!flow.isAllocated()) {
      return;
    }
    auto* mesh = flow.getMesh();
    for (bool lower : {true, false}) {
      const auto side = target(mesh, lower);
      const int y = lower ? mesh->ystart : mesh->yend + 1;
      auto range = lower ? mesh->iterateBndryLowerY() : mesh->iterateBndryUpperY();
      for (; !range.isDone(); ++range) {
        if (*range < mesh->xstart || *range > mesh->xend) {
          continue;
        }
        for (int z = mesh->zstart; z <= mesh->zend; ++z) {
          add(quantity, "present_" + side, 1.0);
          add(quantity, prefix + "_" + side, (lower ? -1.0 : 1.0) * flow(*range, y, z));
        }
      }
    }
  }

  void output(Options& state, const std::string& component, bool diagnose = false) const {
    std::vector<BoutReal> totals;
    for (const auto& [key, value] : local) {
      totals.push_back(value);
    }
    if (totals.empty()) {
      return;
    }
    MPI_Allreduce(MPI_IN_PLACE, totals.data(), static_cast<int>(totals.size()),
                  MPI_DOUBLE, MPI_SUM, BoutComm::get());
    const BoutReal volume = pow(get<BoutReal>(state["rho_s0"]), 3);
    const BoutReal frequency = get<BoutReal>(state["Omega_ci"]);
    const BoutReal density = get<BoutReal>(state["Nnorm"]);
    const BoutReal pressure = SI::qe * density * get<BoutReal>(state["Tnorm"]);
    size_t index = 0;
    for (const auto& [key, value] : local) {
      const auto& [quantity, term] = key;
      const auto total = totals[index++];
      if (term == "rhs" && !diagnose) {
        continue;
      }
      const bool particles = quantity.find("particles_") == 0;
      const bool momentum = quantity.find("momentum_") == 0;
      const bool content = term == "content";
      const bool present = term.find("present_") == 0;
      set_with_attrs(
          state["balance_" + component + "_" + quantity + "_" + term], total,
          {{"time_dimension", "t"},
           {"units", present ? "1"
                             : (particles ? "s^-1"
                                          : (momentum ? (content ? "kg m s^-1" : "N")
                                                      : (content ? "J" : "W")))},
           {"conversion",
            present ? 1.0
                    : (particles
                           ? density
                           : (momentum
                                  ? SI::Mp * density
                                        * std::sqrt(SI::qe * get<BoutReal>(state["Tnorm"])
                                                    / SI::Mp)
                                  : pressure))
                          * volume * (content ? 1.0 : frequency)},
           {"balance_quantity", quantity},
           {"balance_term", term},
           {"source", component},
           {"description", "Global integral over physical cells; simulated domain"}});
    }
  }
};
