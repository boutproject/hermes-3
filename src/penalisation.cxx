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
