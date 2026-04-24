#pragma once
#ifndef EXTERNAL_APAR_H
#define EXTERNAL_APAR_H

#include "component.hxx"

/// Adds an external contribution to the Apar flutter
///
struct ExternalApar : public Component {
  ExternalApar(std::string name, Options& alloptions, Solver* UNUSED(solver));

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
