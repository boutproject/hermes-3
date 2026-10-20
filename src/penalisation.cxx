#include <cmath>
#include <string>

#include <bout/assert.hxx>
#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/field2d.hxx>
#include <bout/field3d.hxx>
#include <bout/globals.hxx>
#include <bout/mesh.hxx>
#include <bout/options.hxx>
#include <bout/output.hxx>
#include <bout/output_bout_types.hxx> // NOLINT
#include <bout/region.hxx>
#include <bout/utils.hxx>

#include "../include/component.hxx"
#include "../include/hermes_utils.hxx"
#include "../include/penalisation.hxx"

using bout::globals::mesh;

Penalisation::Penalisation(const std::string& name, Options& alloptions, Solver*)
    : NamedComponent(name, {
                               readWrite("fields:penalisation_wall"),
                           }) {

  mesh->get(penalisation_wall, "penalisation_wall");
}

void Penalisation::transform_impl(GuardedOptions& state) {

  auto fields = state["fields"];
  set(fields["penalisation_wall"], penalisation_wall);
}

namespace penalisation {
const Field3D neumann(const Field3D& f, const Field3D& ddt_f, const Field3D& pen,
                      const BoutReal pen_time, const BoutReal f_floor) {

  Field3D result{0.0};

  BOUT_FOR(i, f.getRegion("RGN_NOBNDRY")) {
    const BoutReal pen_i = pen[i];
    if (pen_i > 1e-5) {
      const auto ixp = i.xp();
      const auto ixm = i.xm();
      const auto izp = i.zp();
      const auto izm = i.zm();

      // Gradient of pen_A (points toward wall)
      // We want to go OPPOSITE direction (toward plasma)
      const BoutReal dpendx = 0.5 * (pen[ixp] - pen[ixm]); // ∂pen/∂x
      const BoutReal dpendz = 0.5 * (pen[izp] - pen[izm]); // ∂pen/∂z
      // Interior-pointing unit vector = -grad(pen_A)
      const BoutReal grad_mag = sqrt(dpendx * dpendx + dpendz * dpendz);
      BoutReal f_target = 0.0;
      BoutReal weight_sum = 0.0;
      if (grad_mag > 1e-10) {
        // ex, ez: unit vector pointing INTO plasma
        const BoutReal ex = -dpendx / grad_mag;
        const BoutReal ez = -dpendz / grad_mag;

        // Weight each neighbor by how much it aligns
        // with the interior-pointing direction
        // w = max(0, n · e_neighbor)
        const BoutReal w_xp = std::max(0.0, ex);  // +x neighbor
        const BoutReal w_xm = std::max(0.0, -ex); // -x neighbor
        const BoutReal w_zp = std::max(0.0, ez);  // +z neighbor
        const BoutReal w_zm = std::max(0.0, -ez); // -z neighbor

        weight_sum = w_xp + w_xm + w_zp + w_zm;

        if (weight_sum > 1e-10) {
          f_target = (w_xp * f[ixp] + w_xm * f[ixm] + w_zp * f[izp] + w_zm * f[izm])
                     / weight_sum;
        }
      }

      if (weight_sum > 1e-10) {
        result[i] = (1.0 - pen_i) * ddt_f[i] + (pen_i / pen_time) * (f_target - f[i]);
      } else {
        // Deep in wall, no gradient info: damp to floor
        result[i] = (1.0 - pen_i) * ddt_f[i] - (pen_i / pen_time) * (f[i] - f_floor);
      }
    }
  }

  return result;
}
} // namespace penalisation
