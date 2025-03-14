
#ifndef RECALCULATE_METRIC_H
#define RECALCULATE_METRIC_H

#include <bout/bout_types.hxx>

/// Loads geometry and magnetic field quantities from the mesh
///   Rxy, Bpxy, Btxy, hthe, sinty, dx
/// Recalculates co- and contra-variant metric tensors
///
/// Note: Assumes ORTHOGONAL coordinate system.
void recalculate_metric(BoutReal Lnorm, BoutReal Bnorm);

#endif // RECALCULATE_METRIC_H
