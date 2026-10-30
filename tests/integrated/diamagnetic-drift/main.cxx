#include <bout/bout.hxx>
#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/fv_ops.hxx>
#include <bout/mpi_wrapper.hxx>
#include <bout/physicsmodel.hxx>

#include "diamagnetic_drift.hxx"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>

// Evolve only diamagnetic density, pressure and momentum sources. The real
// COMPASS metric and curvature are used without other transport or sources.
class DiamagneticDriftSmoke : public PhysicsModel {
protected:
  int init(bool restarting) override {
    if (restarting or mesh->LocalNz != 1) {
      throw BoutException("Diamagnetic smoke test requires a fresh axisymmetric run");
    }
    SOLVE_FOR3(N, P, M);
    // The grid reader leaves Bxy's processor guards uncommunicated. Match
    // neighbouring cells before the component normalises the curvature.
    mesh->getCoordinates()->setBxy(mesh->getCoordinates()->Bxy());
    volume = mesh->getCoordinates()->cell_volume();
    core_mask = 0.0;
    lower_mask = 0.0;
    upper_mask = 0.0;
    BOUT_FOR_SERIAL(i, N.getRegion("RGN_ALL")) {
      const BoutReal x = static_cast<BoutReal>(mesh->getGlobalXIndex(i.x()) - 2) / 32.0;
      const BoutReal y = TWOPI * mesh->getGlobalYIndexNoBoundaries(i.y()) / 48.0;
      N[i] = 1.0 + 0.2 * x + 0.15 * x * std::sin(y);
      P[i] = 1.5 + 0.2 * x + 0.1 * x * std::cos(y);
      M[i] = 0.3 + 0.05 * x + 0.04 * x * std::sin(2.0 * y);
    }
    if (mesh->firstX() and mesh->periodicY(mesh->xstart)) {
      for (int y = mesh->ystart; y <= mesh->yend; ++y) {
        core_mask(mesh->xstart, y, 0) = 1.0;
      }
    }
    for (auto r = mesh->iterateBndryLowerY(); !r.isDone(); ++r) {
      lower_mask(r.ind, mesh->ystart, 0) = 1.0;
    }
    for (auto r = mesh->iterateBndryUpperY(); !r.isDone(); ++r) {
      upper_mask(r.ind, mesh->yend, 0) = 1.0;
    }
    component_options["units"]["Tesla"] = 1.0;
    component_options["units"]["meters"] = 1.0;
    for (const auto* boundary : {"closed_average", "closed", "open"}) {
      component_options[boundary]["form"] = "divergence";
      component_options[boundary]["radial_boundary"] = boundary;
    }
    averaged = makeComponent("closed_average");
    closed = makeComponent("closed");
    opened = makeComponent("open");
    curvature.covariant = false;
    if (mesh->get(curvature, "bxcv") != 0) {
      throw BoutException("The smoke test requires nonzero grid curvature");
    }
    curvature *= 2.0 / mesh->getCoordinates()->Bxy();

    SAVE_ONCE(volume, core_mask, lower_mask, upper_mask);
    SAVE_REPEAT(Sn, Sp, Sm, reference_error, radial_guard_flux, open_core_flux,
                open_pf_flux, open_outer_flux, unclosed_lower_flux, unclosed_upper_flux);
    return 0;
  }

  int rhs(BoutReal) override {
    mesh->communicate(N, P, M);
    N.applyBoundary();
    P.applyBoundary();
    M.applyBoundary();
    Field3D temperature = P / N;
    // Deliberately unequal target guards exercise the component's T closure.
    for (auto r = mesh->iterateBndryLowerY(); !r.isDone(); ++r) {
      temperature(r.ind, mesh->ystart - 1, 0) *= 3.0;
    }
    for (auto r = mesh->iterateBndryUpperY(); !r.isDone(); ++r) {
      temperature(r.ind, mesh->yend + 1, 0) *= 5.0;
    }
    auto state = makeState(temperature);
    auto closed_state = state.copy();
    auto open_state = state.copy();
    averaged->transform(state);
    closed->transform(closed_state);
    opened->transform(open_state);
    auto& species = state["species"]["d+"];
    Sn = get<Field3D>(species["density_source"]);
    Sp = (2.0 / 3.0) * get<Field3D>(species["energy_source"]);
    Sm = get<Field3D>(species["momentum_source"]);
    ddt(N) = Sn;
    ddt(P) = Sp;
    ddt(M) = Sm;

    // An independent closure: zero the face velocities explicitly, then use
    // the OPEN FV operator. This checks each cell, including each boundary,
    // rather than allowing opposing boundary leaks to cancel in an integral.
    Vector3D velocity = temperature * curvature;
    std::array<BoutReal, 6> probes{};
    auto* coords = mesh->getCoordinates();
    const auto& J = coords->J();
    for (auto r = mesh->iterateBndryLowerY(); !r.isDone(); ++r) {
      const int y = mesh->ystart;
      probes[4] += std::abs(0.125 * (N(r.ind, y, 0) + N(r.ind, y - 1, 0))
                            * (velocity.y(r.ind, y, 0) + velocity.y(r.ind, y - 1, 0))
                            * (J(r.ind, y, 0) + J(r.ind, y - 1, 0))
                            * coords->dx()(r.ind, y, 0) * coords->dz()(r.ind, y, 0));
      velocity.y(r.ind, y - 1, 0) = -velocity.y(r.ind, y, 0);
    }
    for (auto r = mesh->iterateBndryUpperY(); !r.isDone(); ++r) {
      const int y = mesh->yend;
      probes[5] += std::abs(0.125 * (N(r.ind, y, 0) + N(r.ind, y + 1, 0))
                            * (velocity.y(r.ind, y, 0) + velocity.y(r.ind, y + 1, 0))
                            * (J(r.ind, y, 0) + J(r.ind, y + 1, 0))
                            * coords->dx()(r.ind, y, 0) * coords->dz()(r.ind, y, 0));
      velocity.y(r.ind, y + 1, 0) = -velocity.y(r.ind, y, 0);
    }
    for (int y = mesh->ystart; y <= mesh->yend; ++y) {
      if (mesh->firstX()) {
        velocity.x(mesh->xstart - 1, y, 0) = -velocity.x(mesh->xstart, y, 0);
      }
      if (mesh->lastX()) {
        velocity.x(mesh->xend + 1, y, 0) = -velocity.x(mesh->xend, y, 0);
      }
    }
    reference_error = 0.0;
    const std::array<Field3D, 3> quantities{N, P, M};
    const std::array<const char*, 3> names{"density_source", "energy_source",
                                           "momentum_source"};
    for (std::size_t q = 0; q < quantities.size(); ++q) {
      const BoutReal factor = q == 1 ? 2.5 : 1.0;
      Field3D reference = -factor * FV::Div_f_v(quantities[q], velocity, true);
      const auto actual_closed = get<Field3D>(closed_state["species"]["d+"][names[q]]);
      BOUT_FOR_SERIAL(i, N.getRegion("RGN_NOBNDRY")) {
        reference_error[i] =
            std::max(reference_error[i], std::abs(actual_closed[i] - reference[i])
                                             / (1.0 + std::abs(reference[i])));
      }
      // FV deposits radial boundary flux into the physical x guard. Check
      // every such guard for all quantities, with open flux as a control.
      const auto actual_open = get<Field3D>(open_state["species"]["d+"][names[q]]);
      for (int y = mesh->ystart; y <= mesh->yend; ++y) {
        for (const bool inner : {true, false}) {
          if ((inner and !mesh->firstX()) or (!inner and !mesh->lastX())) {
            continue;
          }
          const int x = inner ? mesh->xstart - 1 : mesh->xend + 1;
          probes[0] =
              std::max(probes[0], std::abs(actual_closed(x, y, 0) * volume(x, y, 0)));
          if (q == 0) {
            const std::size_t category =
                !inner ? 3u : (mesh->periodicY(mesh->xstart) ? 1u : 2u);
            probes[category] += std::abs(actual_open(x, y, 0) * volume(x, y, 0));
          }
        }
      }
      if (mesh->firstX() and mesh->periodicY(mesh->xstart)) {
        std::array<BoutReal, 2> local{}, global{};
        for (int y = mesh->ystart; y <= mesh->yend; ++y) {
          local[0] += reference(mesh->xstart, y, 0) * volume(mesh->xstart, y, 0);
          local[1] += volume(mesh->xstart, y, 0);
        }
        bout::globals::mpi->MPI_Allreduce(local.data(), global.data(), 2, MPI_DOUBLE,
                                          MPI_SUM, mesh->getYcomm(mesh->xstart));
        for (int y = mesh->ystart; y <= mesh->yend; ++y) {
          reference(mesh->xstart, y, 0) = global[0] / global[1];
        }
      }
      const auto actual = getNonFinal<Field3D>(species[names[q]]);
      BOUT_FOR_SERIAL(i, N.getRegion("RGN_NOBNDRY")) {
        reference_error[i] =
            std::max(reference_error[i],
                     std::abs(actual[i] - reference[i]) / (1.0 + std::abs(reference[i])));
      }
    }
    std::array<BoutReal, 6> global{};
    bout::globals::mpi->MPI_Allreduce(probes.data(), global.data(), 6, MPI_DOUBLE,
                                      MPI_SUM, BoutComm::get());
    radial_guard_flux = global[0];
    open_core_flux = global[1];
    open_pf_flux = global[2];
    open_outer_flux = global[3];
    unclosed_lower_flux = global[4];
    unclosed_upper_flux = global[5];
    return 0;
  }

private:
  std::unique_ptr<DiamagneticDrift> makeComponent(const char* boundary) {
    auto result =
        std::make_unique<DiamagneticDrift>(boundary, component_options, nullptr);
    result->declareAllSpecies({"d+"});
    return result;
  }
  Options makeState(const Field3D& temperature) {
    Options state;
    auto& species = state["species"]["d+"];
    species["charge"] = 1.0;
    species["temperature"] = temperature;
    species["density"] = N;
    species["pressure"] = P;
    species["momentum"] = M;
    return state;
  }
  Field3D N, P, M, Sn, Sp, Sm, reference_error, core_mask, lower_mask, upper_mask;
  Coordinates::FieldMetric volume;
  VectorMetric curvature;
  Options component_options;
  std::unique_ptr<DiamagneticDrift> averaged, closed, opened;
  BoutReal radial_guard_flux{0.0}, open_core_flux{0.0}, open_pf_flux{0.0},
      open_outer_flux{0.0};
  BoutReal unclosed_lower_flux{0.0}, unclosed_upper_flux{0.0};
};

BOUTMAIN(DiamagneticDriftSmoke);
