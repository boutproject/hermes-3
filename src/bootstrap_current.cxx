
#include "../include/bootstrap_current.hxx"

#include <bout/constants.hxx>
#include <bout/derivs.hxx>
#include <bout/mesh.hxx>
#include <bout/smoothing.hxx>

BootstrapCurrent::BootstrapCurrent(std::string name, Options& alloptions, Solver*) {
  AUTO_TRACE();

  const Options& units = alloptions["units"];

  // Normalisations
  Tnorm = units["eV"];
  Nnorm = units["inv_meters_cubed"];
  Pnorm = SI::qe * Nnorm * Tnorm;
  Bnorm = units["Tesla"];

  Options& options = alloptions[name];

  diagnose =
      options["diagnose"].doc("Save additional diagnostics?").withDefault<bool>(false);

  mesh = bout::globals::mesh;

  // Read Rxy and Btxy from the mesh file
  // Note: These are in SI units
  Field2D Rxy, Btxy;
  GRID_LOAD(Bxy, Rxy, Btxy, J);
  RBt = Rxy * Btxy;

  // This is <1/Bp>, used in fluxSurfaceAverage
  averageJ = averageY(J);
}

void BootstrapCurrent::transform(Options& state) {

  Options& allspecies = state["species"];
  Options& electrons = allspecies["e"];

  // Calculate flux-surface averages
  Field2D pe_av =
      fluxSurfaceAverage(GET_NOBOUNDARY(Field3D, electrons["pressure"])) * Pnorm;

  Field2D Te_av =
      fluxSurfaceAverage(GET_NOBOUNDARY(Field3D, electrons["temperature"])) * Pnorm;

  JparB = fluxSurfaceAverage(get<BoutReal>(electrons["charge"])
                             * GET_NOBOUNDARY(Field3D, electrons["density"])
                             * GET_NOBOUNDARY(Field3D, electrons["velocity"]) * Bxy);

  Field2D pi_av = 0.0;
  Field2D ni_av = 0.0;
  for (auto& kv : allspecies.getChildren()) {
    if (kv.first == "e") {
      continue;
    }

    Options& species = allspecies[kv.first];
    if (!species.isSet("charge")) {
      continue; // Skip uncharged ions
    }
    BoutReal Zi = get<BoutReal>(species["charge"]);
    if (fabs(Zi) < 1e-3) {
      continue;
    }

    pi_av += fluxSurfaceAverage(GET_NOBOUNDARY(Field3D, species["pressure"])) * Pnorm;
    ni_av += fluxSurfaceAverage(GET_NOBOUNDARY(Field3D, species["density"])) * Nnorm;

    JparB += fluxSurfaceAverage(Zi * GET_NOBOUNDARY(Field3D, species["density"])
                                * GET_NOBOUNDARY(Field3D, species["velocity"]) * Bxy);
  }

  Field2D p_av = pe_av + pi_av;
  Field2D Ti_av = pi_av / ni_av;

  // Communicate and set boundaries
  mesh->communicate(p_av, pe_av, Ti_av, Te_av);

  // Radial gradients
  Field2D dp_dpsi = DDX(p_av);
  Field2D dlnTe_dpsi = DDX(log(Te_av));
  Field2D dlnTi_dpsi = DDX(log(Ti_av));

  // Electron collisionality
  // nu_estar = 0.012 n_20 Zeff qR / epsilon^3/2 T_ekev^2
  Field2D nu_estar = 0.01;
  Field2D nu_istar = 0.1;

  // Trapped fraction f_t
  Field2D f_t = 0.5;

  // Effective ion charge
  Field2D Zeff = 1.0;

  // Loop over flux surfaces
  for (int ix = mesh->xstart; ix <= mesh->xend; ++ix) {
    if (!mesh->periodicY(ix)) {
      // Not periodic => Outside core
      break;
    }
    // Periodic Y i.e. in the core

    const BoutReal f = f_t(ix, mesh->ystart);
    const BoutReal Z = Zeff(ix, mesh->ystart);
    const BoutReal nu_e = nu_estar(ix, mesh->ystart);
    const BoutReal nu_i = nu_istar(ix, mesh->ystart);

    // Calculate coefficients
    const BoutReal L31_ix = L31(Z, f, nu_e);
    const BoutReal L32_ix = L32(Z, f, nu_e);
    const BoutReal L34_ix = L34(Z, f, nu_e);
    const BoutReal alpha_ix = alpha(f, nu_i);

    const BoutReal pe_ix = pe_av(ix, mesh->ystart);
    const BoutReal dp_dpsi_ix = dp_dpsi(ix, mesh->ystart);
    const BoutReal dlnTe_dpsi_ix = dlnTe_dpsi(ix, mesh->ystart);
    const BoutReal dlnTi_dpsi_ix = dlnTi_dpsi(ix, mesh->ystart);

    const BoutReal RBt_ix = RBt(ix, mesh->ystart);

    // Steady state parallel current <J||B>
    BoutReal JparB_ix =
        -RBt_ix
        * (L31_ix * dp_dpsi_ix
           + pe_ix * (L32_ix * dlnTe_dpsi_ix + L34_ix * alpha_ix * dlnTi_dpsi_ix));

    // Electron momentum source term that drives current toward
    // steady state solution. Note: Should have zero parallel divergence
    // so Force / B = const on flux surface

    for (int jy = mesh->ystart; jy <= mesh->yend; ++jy) {
      JparB_bs(ix, jy) = JparB_ix;
    }
  }
}

void BootstrapCurrent::outputVars(Options& state) {
  if (!diagnose) {
    return;
  }

  const auto Nnorm = get<BoutReal>(state["Nnorm"]);
  const auto Tnorm = get<BoutReal>(state["Tnorm"]);
  const BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation
  const auto Bnorm = get<BoutReal>(state["Bnorm"]);
  const auto Cs0 = get<BoutReal>(state["Cs0"]);

  set_with_attrs(state["JparB"], JparB,
                 {{"time_dimension", "t"},
                  {"units", "s-1"},
                  {"conversion", SI::qe * Cs0 * Nnorm * Bnorm},
                  {"standard_name", "<J||B>"},
                  {"long_name", "Flux surface average <J dot B>"},
                  {"source", "bootstrap_current"}});
}

// Equation 14a
BoutReal BootstrapCurrent::F31(BoutReal X, BoutReal Z) {
  const BoutReal X2 = SQ(X);
  const BoutReal X3 = X2 * X;
  const BoutReal X4 = SQ(X2);
  return (1. + 1.4 / (Z + 1.)) * X - 1.9 / (Z + 1.) * X2 + 0.3 / (Z + 1.) * X3
         + 0.2 / (Z + 1.) * X4;
}

BoutReal BootstrapCurrent::L31(BoutReal Z, BoutReal f, BoutReal nu_e) {
  const BoutReal f31_teff =
      f / (1. + (1. - 0.1 * f) * sqrt(nu_e) + 0.5 * (1. - f) * nu_e / Z);
  return F31(f31_teff, Z);
}

BoutReal BootstrapCurrent::L32(BoutReal Z, BoutReal f, BoutReal nu_e) {
  // Equation 15d
  const BoutReal f32ee_teff =
      f / (1. + 0.26 * (1. - f) * sqrt(nu_e) + 0.18 * (1. - 0.37 * f) * nu_e / sqrt(Z));

  // Equation 15b
  const BoutReal X = f32ee_teff;
  const BoutReal X2 = SQ(f32ee_teff);
  const BoutReal X3 = X2 * X;
  const BoutReal X4 = SQ(X2);
  const BoutReal F32_ee = (0.05 + 0.62 * Z) / (Z * (1. + 0.44 * Z)) * (X - X4)
                          + 1. / (1. + 0.22 * Z) * (X2 - X4 - 1.2 * (X3 - X4))
                          + 1.2 / (1. + 0.5 * Z) * X4;

  // Equation 15e
  const BoutReal f32ei_teff =
      f / (1. + (1. + 0.6 * f) * sqrt(nu_e) + 0.85 * (1. - 0.37 * f) * nu_e * (1. + Z));

  // Equation 15c
  const BoutReal Y = f32ei_teff;
  const BoutReal Y2 = SQ(Y);
  const BoutReal Y3 = Y2 * Y;
  const BoutReal Y4 = SQ(Y2);
  const BoutReal F32_ei = -(0.56 + 1.93 * Z) / (Z * (1. + 0.44 * Z)) * (Y - Y4)
                          + 4.95 / (1. + 2.48 * Z) * (Y2 - Y4 - 0.55 * (Y3 - Y4))
                          - 1.2 / (1. + 0.5 * Z) * Y4;

  // Equation 15a
  return F32_ee + F32_ei;
}

BoutReal BootstrapCurrent::L34(BoutReal Z, BoutReal f, BoutReal nu_e) {
  // Equation 16b
  const BoutReal f34_teff =
      f / (1. + (1. - 0.1 * f) * sqrt(nu_e) + 0.5 * (1. - 0.5 * f) * nu_e / Z);

  // Equation 16a
  return F31(f34_teff, Z);
}

BoutReal BootstrapCurrent::alpha(BoutReal f, BoutReal nu_i) {
  // Equation 17a
  const BoutReal f2 = SQ(f);
  const BoutReal alpha_0 = -1.17 * (1. - f) / (1. - 0.22 * f - 0.19 * f2);

  // Equation 17b
  // Note: The sign of the 0.315 coefficient is corrected in the
  // O.Sauter Erratum DOI: 10.1063/1.1517052
  const BoutReal f6 = SQ(f2 * f);
  return ((alpha_0 + 0.25 * (1. - f2) * sqrt(nu_i)) / (1. + 0.5 * sqrt(nu_i))
          + 0.315 * SQ(nu_i) * f6)
         / (1. + 0.15 * SQ(nu_i) * f6);
}

Field2D BootstrapCurrent::fluxSurfaceAverage(const Field3D& f) {
  return fluxSurfaceAverage(DC(f));
}

Field2D BootstrapCurrent::fluxSurfaceAverage(const Field2D& f) {
  return averageY(f * J) / averageJ;
}

Field2D BootstrapCurrent::trappedFraction() {
  // Maximum magnetic field on each flux surface
  const Field2D Bmax = maxY(Bxy);

  // Integrate λ / < sqrt(1 - λ) >
  // over 0 < λ < 1 / Bmax
  //
  // Transform x = 2λBmax - 1
  // and integrate:
  //
  // (1/(Bmax^2 sqrt(2))) (1 + x)/sqrt(1 - x) sqrt(1 - x) / <sqrt(1 - x B/Bmax)>
  //                      |------ w(x) -----|
  // over the range -1 < x < 1
  //
  // The sqrt(1 - x) factor removes the singularity, moving it into
  // the weight function w(x)
  //
  // Use Gauss-Jacobi rule with 5 points. alpha = -0.5, beta = 1
  // Python:
  //     from scipy.special import roots_jacobi
  //     x, w = roots_jacobi(5, -0.5, 1.0)
  //
  const std::vector<std::pair<BoutReal, BoutReal>> quadrature = {
      {-0.7856692692946651, 0.05558088292849234},
      {-0.3424372137469275, 0.2938121661741237},
      {0.19893554984718564, 0.7206376843619846},
      {0.6807500544226858, 1.194533824278887},
      {0.9627065930574353, 1.5066716085847658}};

  Field2D integral = 0.0;
  for (const auto& [x, w] : quadrature) {
    integral += w * sqrt(1. - x) / fluxSurfaceAverage(sqrt(1. - x * Bxy / Bmax));
  }

  // Trapped fraction
  return 1. - (3. / 4) * fluxSurfaceAverage(SQ(Bxy)) / (SQ(Bmax) * sqrt(2.)) * integral;
}

Field2D BootstrapCurrent::maxY(const Field2D& f) {
  TRACE("maxY(Field2D)");

  Mesh* mesh = f.getMesh();
  int ngx = mesh->LocalNx;
  int ngy = mesh->LocalNy;

  Array<BoutReal> input(ngx), result(ngx);

  // Maximum on this processor
  for (int x = 0; x < ngx; x++) {
    input[x] = f(x, mesh->ystart);
    // Maximum value, not including boundaries
    for (int y = mesh->ystart; y <= mesh->yend; y++) {
      input[x] = std::max(input[x], f(x, y));
    }
  }

  Field2D r{emptyFrom(f)};

  /// NOTE: This only works if there are no branch-cuts
  MPI_Comm comm_inner = mesh->getYcomm(0);

  int np;
  MPI_Comm_size(comm_inner, &np);

  if (np == 1) {
    for (int x = 0; x < ngx; x++) {
      for (int y = 0; y < ngy; y++) {
        r(x, y) = input[x];
      }
    }
  } else {
    MPI_Allreduce(input.begin(), result.begin(), ngx, MPI_DOUBLE, MPI_MAX, comm_inner);
    for (int x = 0; x < ngx; x++) {
      for (int y = 0; y < ngy; y++) {
        r(x, y) = result[x];
      }
    }
  }
  return r;
}
