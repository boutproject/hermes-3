#pragma once
#ifndef FREE_STREAMING_H
#define FREE_STREAMING_H

#include <bout/coordinates.hxx>
#include <bout/difops.hxx>
#include <bout/field2d.hxx>
#include <bout/field3d.hxx>
#include <bout/globals.hxx>
#include <bout/mesh.hxx>
#include <bout/smoothing.hxx>

/// Field line length along y (normalised), from the metric.
/// shortcut: averageY needs a grid without branch cuts (e.g. closed field lines);
/// on other grids give the length as an input instead
inline Field2D fieldLineLength() {
  Mesh* mesh = bout::globals::mesh;
  const Coordinates* coord = mesh->getCoordinates();
  return averageY(Field2D{sqrt(coord->g_22()) * coord->dy()})
         * mesh->GlobalNyNoBoundaries;
}

/// 1 / (1 + lambda / (c L_par)) with 1/L_par = 1/L_line, plus |Grad_par T| / T
/// if \p gradient.
///
/// Reduces a Braginskii coefficient, which assumes lambda = v_th tau << L_par, smoothly
/// when the plasma is collisionless. L_par is the field line length where T is flat
/// along B, and the local gradient length where T is steep.
/// With \p gradient, Grad_par is only set inside the domain: communicate the result
/// before using guard cells.
inline Field3D freeStreamingFactor(const Field3D& T, const Field3D& lambda,
                                   const Field2D& L_line, BoutReal c, bool gradient) {
  Field3D inv_L_par{1. / L_line};
  if (gradient) {
    inv_L_par += abs(Grad_par(T)) / floor(T, 1e-5);
  }
  // floor: Newton trial states can have tau < 0 (n < 0), where 1 + lambda / (c L_par)
  // would cross zero; keeps the factor in (0, 1] and leaves the coefficient unchanged there
  return 1. / (1. + floor(lambda, 0.0) * inv_L_par / c);
}

#endif // FREE_STREAMING_H
