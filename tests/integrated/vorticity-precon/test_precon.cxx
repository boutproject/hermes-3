#include <bout/petsclib.hxx>

#include "component_scheduler.hxx"
#include "electromagnetic.hxx"
#include "evolve_momentum.hxx"
#include "fixed_density.hxx"
#include "vorticity.hxx"

#include <bout/bout.hxx>
#include <bout/constants.hxx>
#include <bout/derivs.hxx>
#include <bout/difops.hxx>
#include <bout/field_factory.hxx>
#include <bout/fv_ops.hxx>
#include <bout/solver.hxx>

#include <iostream>
#include <limits>

namespace {
class ResidualSolver : public Solver {
public:
  int run() override { return 0; }
  Field3D& residual(const std::string& name) {
    for (const auto& field : f3d) {
      if (field.name == name) {
        return *field.F_var;
      }
    }
    throw BoutException("Missing test field {}", name);
  }
};

void checkBlock(bool momentum_first, bool em) {
  auto* mesh = bout::globals::mesh;
  auto& options = Options::root();
  Options schedule;
  std::string order = momentum_first ? "e,vorticity" : "vorticity,e";
  schedule["components"] = em ? "i," + order + ",em" : order;
  schedule["autosort"] = false;
  ResidualSolver solver;
  ComponentScheduler components(schedule, options, &solver);

  const BoutReal N = em ? 3e-6 : 1.25;
  const BoutReal Ni = 2e-6;
  const BoutReal Z = -1.0;
  const BoutReal AA = options["e"]["AA"];
  Options state;
  state["time"] = 0.;
  state["species"]["e"] = {
      {"AA", AA}, {"charge", Z}, {"density", Field3D{N}}, {"temperature", Field3D{1.}}};
  if (em) {
    state["species"]["i"] = {{"momentum", Field3D{0.}}, {"velocity", Field3D{0.}}};
  }
  components.transform(state);

  auto make = [&](const std::string& expr) {
    return FieldFactory::get()->create3D(expr, &options, mesh);
  };
  Field3D phi = make("sin(2*pi*x + 0.2)*cos(y+0.3)");
  mesh->communicate(phi);
  phi.setBoundary("test_phi");
  phi.applyBoundary();
  const Field3D vort = Delp2(phi);
  Field3D momentum = make("0.3*cos(3*x+y)");
  Field3D current;
  if (em) {
    // Manufacture the canonical momentum from the configured Ampère operator.
    // On this mesh its nested solver drops G1*d/dx and keeps g11*d2/dx2.
    Field3D apar = make("cos(2*x+0.1)*sin(2*y+0.3)");
    mesh->communicate(apar);
    apar.applyBoundary("dirichlet");
    const BoutReal beta =
        SI::mu0 * SI::qe * options["units"]["inv_meters_cubed"].as<BoutReal>();
    const BoutReal alpha = 1e-5 * (1. / AA + 1.);
    momentum =
        (AA / Z) * (alpha * apar - mesh->getCoordinates()->g11() * D2DX2(apar) / beta);
    current = (Z / AA) * (momentum - Z * N * apar) - Ni * apar;
  } else {
    current = (Z / AA) * momentum;
  }
  mesh->communicate(current);
  current.applyBoundary("dirichlet"); // Fixed current at the target face

  const Field3D ac = make("cos(z)");
  // Use exactly zero toroidal mean independently of the configured z extent.
  const Field3D ac_part = 0.01 * (ac - DC(ac));
  Field3D& rV = solver.residual("Vort");
  Field3D& rP = solver.residual("NVe");
  auto check = [&](bool viscosity, BoutReal gamma = 0.7) {
    Options callback_state = state.copy();
    rV = vort - gamma * Div_par(current) + ac_part;
    rP = momentum + gamma * Z * N * Grad_par(phi) + ac_part;
    if (viscosity) {
      const Field3D eta{0.25};
      callback_state["species"]["e"]["eta_par"] = eta;
      callback_state["species"]["e"]["sheath_velocity_closure"] = Field3D{-1.0};
      Field3D velocity = (momentum + ac_part) / (AA * N);
      mesh->communicate(velocity);
      velocity.applyBoundary("dirichlet");
      rP -= gamma * FV::Div_par_K_Grad_par(eta, velocity);
    }

    // The integrator supplies interior momentum residuals and may evolve Vort boundaries.
    for (auto i : rP.getRegion("RGN_BNDRY")) {
      rP[i] = std::numeric_limits<BoutReal>::quiet_NaN();
      rV[i] = 7.0;
    }
    components.precon(callback_state, gamma);
    const BoutReal v_error = max(abs(rV - vort - ac_part), true, "RGN_NOBNDRY");
    const BoutReal p_error = max(abs(rP - momentum - ac_part), true, "RGN_NOBNDRY");
    const BoutReal b_error = max(abs(rV - 7.0), true, "RGN_BNDRY");
    if (!(v_error < 1e-10 && p_error < 1e-11 && b_error == 0.0)) {
      throw BoutException("n=0 block (momentum_first={}, em={}, viscosity={}): Vort {}, "
                          "NV {}, boundary {}",
                          momentum_first, em, viscosity, v_error, p_error, b_error);
    }
  };
  check(false);
  if (!em) {
    // Keep gamma and time fixed across callbacks so a viscosity switch must
    // invalidate the cached factorisation even within a single timestep.
    check(true);
    check(false);
    check(true);
    check(true, 0.0);
  }
}
} // namespace

int main(int argc, char** argv) {
  if (BoutInitialise(argc, argv) != 0) {
    return 1;
  }
  int status = 0;
  try {
    // Hermes normally installs these as typed values before creating components.
    Options::root()["units"] = {{"meters", 1.},
                                {"seconds", 1.},
                                {"Tesla", 1.},
                                {"eV", 1.},
                                {"inv_meters_cubed", 5e29}};
    PetscLib petsc; // Keep PETSc alive across independent scheduler instances.
    for (bool em : {false, true}) {
      checkBlock(true, em);
      checkBlock(false, em);
    }
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    status = 1;
  }
  BoutFinalise(false);
  return status;
}
