#pragma once
#ifndef REACTION_DIAGNOSTIC_H
#define REACTION_DIAGNOSTIC_H

#include <functional>
#include <string>

#include <bout/bout_enum_class.hxx>
#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/constants.hxx>

#include "component.hxx"

BOUT_ENUM_CLASS(ReactionDiagnosticType, density, momentum, energy);
std::string toString(ReactionDiagnosticType diag_type);

typedef std::function<Field3D(const Field3D&)> DiagnosticTransformerType;

// Some pre-defined functions to use as diagnostic transformers. Default is identity.
inline Field3D identity(const Field3D& src_fld) { return src_fld; };
inline Field3D negate(const Field3D& src_fld) { return -src_fld; };

struct ReactionDiagnostic {
  static inline std::string default_std_name(ReactionDiagnosticType type) {
    std::string standard_name;
    switch (type) {
    case ReactionDiagnosticType::density:
      standard_name = "particle source";
      break;
    case ReactionDiagnosticType::momentum: // fall through
    case ReactionDiagnosticType::energy:
      standard_name = fmt::format("{:s} transfer", toString(type));
      break;
    default:
      throw BoutException("ReactionDiagnostic type not set up!");
    }
    return standard_name;
  }
  /**
   * @todo just use a default param for standard_name?
   *
   */
  ReactionDiagnostic(const std::string& name, const std::string& long_name,
                     ReactionDiagnosticType type, const std::string& source,
                     DiagnosticTransformerType transformer = identity)
      : ReactionDiagnostic(name, long_name, type, source, default_std_name(type),
                           transformer){};

  ReactionDiagnostic(const std::string& name, const std::string& long_name,
                     ReactionDiagnosticType type, const std::string& source,
                     const std::string& standard_name,
                     DiagnosticTransformerType transformer)
      : long_name(long_name), name(name), source(source), standard_name(standard_name),
        transformer(transformer), type(type) {

    // Extract normalisation units from the root options
    Options& options = Options::root();
    const BoutReal Nnorm = get<BoutReal>(options["units"]["inv_meters_cubed"]);
    const BoutReal Tnorm = get<BoutReal>(options["units"]["eV"]);
    const BoutReal Omega_ci = 1 / get<BoutReal>(options["units"]["seconds"]);

    switch (type) {
    case ReactionDiagnosticType::density: {
      this->conversion = Nnorm * Omega_ci;
      this->units = "m^-3 s^-1";
      break;
    }
    case ReactionDiagnosticType::momentum: {
      const BoutReal Cs0 = std::sqrt(SI::qe * Tnorm / SI::Mp);
      this->conversion = SI::Mp * Nnorm * Cs0 * Omega_ci;
      this->units = "kg m^-2 s^-2";
      break;
    }
    case ReactionDiagnosticType::energy: {
      const BoutReal Pnorm = SI::qe * Tnorm * Nnorm; // Pressure normalisation
      this->conversion = Pnorm * Omega_ci;
      this->units = "W / m^3";
      break;
    }
    default:
      throw BoutException("ReactionDiagnostic type not set up!");
    }
  }

  void set_attrs(Options& state) {
    state[this->name].setAttributes({{"time_dimension", "t"},
                                     {"units", this->units},
                                     {"conversion", this->conversion},
                                     {"standard_name", this->standard_name},
                                     {"long_name", this->long_name},
                                     {"source", this->source}});
  }

  Field3D transform(Field3D src_fld) { return this->transformer(src_fld); }

  const std::string name;
  const std::string long_name;
  const std::string source;
  const std::string standard_name;
  const ReactionDiagnosticType type;

private:
  std::string conversion;
  const DiagnosticTransformerType transformer;
  std::string units;
};

#endif