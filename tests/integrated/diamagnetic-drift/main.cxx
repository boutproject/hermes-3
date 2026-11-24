#include <bout/bout.hxx>
#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/mpi_wrapper.hxx>
#include <bout/physicsmodel.hxx>

#include "diamagnetic_drift.hxx"

#include <array>
#include <cmath>
#include <memory>

// Evolve only production diamagnetic transport on the COMPASS grid.
class DiamagneticDriftSmoke : public PhysicsModel {
protected:
  int init(bool restarting) override {
    if (restarting) {
      throw BoutException("Diamagnetic smoke test requires a fresh run");
    }
    SOLVE_FOR3(N, P, M);
    mesh->getCoordinates()->setBxy(mesh->getCoordinates()->Bxy());
    volume = mesh->getCoordinates()->cell_volume();
    core_mask = 0.0;
    BOUT_FOR_SERIAL(i, N.getRegion("RGN_ALL")) {
      const BoutReal x = static_cast<BoutReal>(mesh->getGlobalXIndex(i.x()) - 2) / 32.0;
      const BoutReal y = TWOPI * mesh->getGlobalYIndexNoBoundaries(i.y()) / 48.0;
      const BoutReal z = TWOPI * i.z() / mesh->LocalNz;
      N[i] = 1.0 + 0.2 * x + 0.15 * x * std::sin(y) + 0.02 * x * std::cos(z);
      P[i] = 1.5 + 0.2 * x + 0.1 * x * std::cos(y);
      M[i] = 0.3 + 0.05 * x + 0.04 * x * std::sin(2.0 * y);
      if (mesh->firstX() and i.x() == mesh->xstart and mesh->periodicY(i.x())
          and i.y() >= mesh->ystart and i.y() <= mesh->yend) {
        core_mask[i] = 1.0;
      }
    }
    component_options["units"]["Tesla"] = 1.0;
    component_options["units"]["meters"] = 1.0;
    averaged = makeComponent("closed_average");
    closed = makeComponent("closed");
    opened = makeComponent("open");
    SAVE_ONCE(volume, core_mask);
    SAVE_REPEAT(Sn, Sp, Sm, Sn_closed, Sp_closed, Sm_closed, radial_guard_flux,
                open_radial_flux);
    return 0;
  }

  int rhs(BoutReal) override {
    mesh->communicate(N, P, M);
    N.applyBoundary();
    P.applyBoundary();
    M.applyBoundary();
    Field3D temperature = P / N;
    for (auto r = mesh->iterateBndryLowerY(); !r.isDone(); ++r) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        temperature(r.ind, mesh->ystart - 1, z) *= 3.0;
      }
    }
    for (auto r = mesh->iterateBndryUpperY(); !r.isDone(); ++r) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        temperature(r.ind, mesh->yend + 1, z) *= 5.0;
      }
    }
    auto state = makeState(temperature);
    auto closed_state = state.copy(), open_state = state.copy();
    averaged->transform(state);
    closed->transform(closed_state);
    opened->transform(open_state);
    auto& species = state["species"]["d+"];
    auto& closed_species = closed_state["species"]["d+"];
    Sn = get<Field3D>(species["density_source"]);
    Sp = (2.0 / 3.0) * get<Field3D>(species["energy_source"]);
    Sm = get<Field3D>(species["momentum_source"]);
    Sn_closed = get<Field3D>(closed_species["density_source"]);
    Sp_closed = (2.0 / 3.0) * get<Field3D>(closed_species["energy_source"]);
    Sm_closed = get<Field3D>(closed_species["momentum_source"]);
    ddt(N) = Sn;
    ddt(P) = Sp;
    ddt(M) = Sm;
    // Check production radial deposition directly; no copied flux operator.
    std::array<BoutReal, 2> local{}, global{};
    for (const auto* name : {"density_source", "energy_source", "momentum_source"}) {
      const auto closed_source = get<Field3D>(closed_species[name]);
      const auto open_source = get<Field3D>(open_state["species"]["d+"][name]);
      for (bool inner : {true, false}) {
        if ((inner and !mesh->firstX()) or (!inner and !mesh->lastX())) {
          continue;
        }
        const int x = inner ? mesh->xstart - 1 : mesh->xend + 1;
        for (int y = mesh->ystart; y <= mesh->yend; ++y) {
          for (int z = 0; z < mesh->LocalNz; ++z) {
            local[0] += std::abs(closed_source(x, y, z) * volume(x, y, z));
            local[1] += std::abs(open_source(x, y, z) * volume(x, y, z));
          }
        }
      }
    }
    bout::globals::mpi->MPI_Allreduce(local.data(), global.data(), 2, MPI_DOUBLE, MPI_SUM,
                                      BoutComm::get());
    radial_guard_flux = global[0];
    open_radial_flux = global[1];
    return 0;
  }

private:
  std::unique_ptr<DiamagneticDrift> makeComponent(const char* boundary) {
    component_options[boundary]["form"] = "divergence";
    component_options[boundary]["radial_boundary"] = boundary;
    auto result =
        std::make_unique<DiamagneticDrift>(boundary, component_options, nullptr);
    result->declareAllSpecies({"d+"});
    return result;
  }
  Options makeState(const Field3D& temperature) {
    return {{"species",
             {{"d+",
               {{"charge", 1.0},
                {"temperature", temperature},
                {"density", N},
                {"pressure", P},
                {"momentum", M}}}}}};
  }
  Field3D N, P, M, Sn, Sp, Sm, Sn_closed, Sp_closed, Sm_closed, core_mask;
  Coordinates::FieldMetric volume;
  Options component_options;
  std::unique_ptr<DiamagneticDrift> averaged, closed, opened;
  BoutReal radial_guard_flux{0.0}, open_radial_flux{0.0};
};

BOUTMAIN(DiamagneticDriftSmoke);
