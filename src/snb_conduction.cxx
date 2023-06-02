#include <bout/constants.hxx>
#include "../include/snb_conduction.hxx"

#include <bout/bout.hxx>
using bout::globals::mesh;

void SNBConduction::transform(Options& state) {
  auto units = state["units"];
  const auto rho_s0 = get<BoutReal>(units["meters"]);
  const auto Tnorm = get<BoutReal>(units["eV"]);
  const auto Nnorm = get<BoutReal>(units["inv_meters_cubed"]);
  const auto Omega_ci = 1. / get<BoutReal>(units["seconds"]);

  Options& electrons = state["species"]["e"];
  // Note: Needs boundary conditions on temperature
  const Field3D Te = GET_VALUE(Field3D, electrons["temperature"]) * Tnorm; // eV
  const Field3D Ne = GET_VALUE(Field3D, electrons["density"]) * Nnorm;     // In m^-3

  // SNB non-local heat flux. Also returns the Spitzer-Harm value for comparison
  // Note: Te in eV, Ne in Nnorm
  Field3D dy_orig = mesh->getCoordinates()->dy;
  mesh->getCoordinates()->dy *= rho_s0; // Convert distances to m

  // Inputs in eV and m^-3
  Div_Q_SNB = snb.divHeatFlux(Te, Ne, &Div_Q_SH);

  // Restore the metric tensor
  mesh->getCoordinates()->dy = dy_orig;

  // Normalise from eV/m^3/s
  Div_Q_SNB /= Tnorm * Nnorm * Omega_ci;
  Div_Q_SH /= Tnorm * Nnorm * Omega_ci;

  // Divergence of heat flux appears with a minus sign in energy source:
  //
  // ddt(3/2 P) = .. - Div(q)
  subtract(electrons["energy_source"], Div_Q_SNB);
}

void SNBConduction::outputVars(Options& state) {
  AUTO_TRACE();

  if (diagnose) {
    auto Nnorm = get<BoutReal>(state["Nnorm"]);
    auto Tnorm = get<BoutReal>(state["Tnorm"]);
    auto Omega_ci = get<BoutReal>(state["Omega_ci"]);

    BoutReal DivQnorm = SI::qe * Tnorm * Nnorm * Omega_ci;

    set_with_attrs(state["Div_Q_SH"], Div_Q_SH,
                   {{"time_dimension", "t"},
                    {"units", "W / m^3"},
                    {"conversion", DivQnorm},
                    {"long_name", "Divergence of Spitzer-Harm electron heat conduction"},
                    {"source", "snb_conduction"}});

    set_with_attrs(state["Div_Q_SNB"], Div_Q_SNB,
                   {{"time_dimension", "t"},
                    {"units", "W / m^3"},
                    {"conversion", DivQnorm},
                    {"long_name", "Divergence of SNB electron heat conduction"},
                    {"source", "snb_conduction"}});
  }
}
