#pragma once
#ifndef EXTERNAL_APAR_H
#define EXTERNAL_APAR_H

#include "component.hxx"
#include "guarded_options.hxx"

#include <bout/field3d.hxx>
#include <bout/options.hxx>

#include <string>

/// Adds an external contribution to the Apar flutter
///
struct ExternalApar : public Component {
  ExternalApar(std::string name, Options& alloptions, [[maybe_unused]] Solver* solver);

  /// Saves the added field to output
  void outputVars(Options& state) override;

private:
  /// Adds to the Apar_flutter field
  ///
  /// - fields
  ///   - Apar_flutter
  void transform_impl(GuardedOptions& state) override;

  Field3D external_apar; ///< The external field
};

namespace {
RegisterComponent<ExternalApar> registercomponentexternalapar("external_apar");
}

#endif // EXTERNAL_APAR_H
