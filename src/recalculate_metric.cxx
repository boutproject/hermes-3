#include <bout/field2d.hxx>
#include <bout/globals.hxx>
#include <bout/mesh.hxx>
#include <bout/output.hxx>
#include <bout/utils.hxx>

#include "../include/recalculate_metric.hxx"

using bout::globals::mesh;

void recalculate_metric(BoutReal Lnorm, BoutReal Bnorm) {
  // Load Rxy etc from from the mesh
  Field2D Rxy, Bpxy, Btxy, hthe, sinty;
  GRID_LOAD5(Rxy, Bpxy, Btxy, hthe, sinty); // Load metrics

  Coordinates* coord = mesh->getCoordinates();

  // Checking for dpsi and qinty used in BOUT grids
  Field2D dx;
  if (!mesh->get(dx, "dpsi")) {
    output << "\tUsing dpsi as the x grid spacing\n";
    coord->setDx(dx); // Only use dpsi if found
  } else {
    // dx will have been read already from the grid
    output << "\tUsing dx as the x grid spacing\n";
  }

  Rxy /= Lnorm;
  hthe /= Lnorm;
  sinty *= SQ(Lnorm) * Bnorm;
  coord->setDx(coord->dx() / (SQ(Lnorm) * Bnorm));

  Bpxy /= Bnorm;
  Btxy /= Bnorm;
  coord->setBxy(coord->Bxy() / (Bnorm));

  // Calculate metric components
  if (Options::root()["mesh"]["paralleltransform"]["type"].as<std::string>()
      == "shifted") {
    sinty = 0.0; // I disappears from metric
  }

  BoutReal sbp = 1.0; // Sign of Bp
  if (min(Bpxy, true) < 0.0) {
    sbp = -1.0;
  }

  const Field2D g11 = SQ(Rxy * Bpxy);
  const Field2D g22 = 1.0 / SQ(hthe);
  const auto g33 = SQ(sinty) * g11 + SQ(coord->Bxy()) / g11;
  const Field2D g12 = 0.0;
  const Field2D g13 = -sinty * g11;
  const Field2D g23 = -sbp * Btxy / (hthe * Bpxy * Rxy);

  const Field2D g_11 = 1.0 / g11 + SQ(sinty * Rxy);
  const auto g_22 = SQ(coord->Bxy() * hthe / Bpxy);
  const Field2D g_33 = Rxy * Rxy;
  const Field2D g_12 = sbp * Btxy * hthe * sinty * Rxy / Bpxy;
  const Field2D g_13 = sinty * Rxy * Rxy;
  const Field2D g_23 = sbp * Btxy * hthe * Rxy / Bpxy;

  coord->setMetricTensor(ContravariantMetricTensor(g11, g22, g33, g12, g13, g23),
                         CovariantMetricTensor(g_11, g_22, g_33, g_12, g_13, g_23));

  const Field2D J = hthe / Bpxy;
  coord->setJ(J);
}
