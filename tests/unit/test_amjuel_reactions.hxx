#ifndef TEST_AMJUEL_REACTIONS_H__
#define TEST_AMJUEL_REACTIONS_H__

#include "test_reactions.hxx"

#include "../../include/amjuel_helium.hxx"
#include "../../include/amjuel_hyd_ionisation.hxx"
#include "../../include/amjuel_hyd_recombination.hxx"

/**
 * @brief Base fixture for tests of AmjuelReaction subclasses.
 *
 * @tparam RTYPE The reaction class; must derive from AmjuelReaction.
 */
template <typename RTYPE>
class AmjuelReactionTest : public ReactionTest<RTYPE> {

  static_assert(std::is_base_of<AmjuelReaction, RTYPE>(),
                "Template arg to AmjuelReactionTest must derive from AmjuelReaction");

public:
  AmjuelReactionTest(std::string lbl, std::string reaction_str,
                     std::string heavy_reactant)
      : ReactionTest<RTYPE>(lbl, reaction_str), heavy_reactant(heavy_reactant) {}

protected:
  /// Name of the ion/atom reactant
  std::string heavy_reactant;

  /**
   * Set up test state for tests of Amjuel-based reactions
   *
   */
  Options generate_state() override {
    // Clunky way to record the ion and neutral species names
    std::string atom(heavy_reactant);
    atom.erase(std::remove(atom.begin(), atom.end(), '+'), atom.end());
    std::string ion = atom + "+";

    const std::map<std::string, BoutReal> heavy_reactant_masses = {
        {"h", 1.0}, {"he", 4.0}, {"d", 2.0}, {"t", 3.0}, {"e", 1. / 1836}};
    const BoutReal atom_charge = 0.0;
    const BoutReal ion_charge = 1.0;
    std::string comp_name("test" + this->lbl);
    Options state{{comp_name, {{"type", this->reaction_str}}},
                  {"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"species",
                   {{"e", {{"AA", heavy_reactant_masses.at("e")}, {"velocity", 1.0}}},
                    {atom,
                     {{"AA", heavy_reactant_masses.at(atom)},
                      {"charge", atom_charge},
                      {"density", 1.0},
                      {"temperature", 1.0},
                      {"velocity", 1.0}}},
                    {ion,
                     {{"AA", heavy_reactant_masses.at(atom)},
                      {"charge", ion_charge},
                      {"density", 1.0},
                      {"temperature", 1.0},
                      {"velocity", 1.0}}}}}};

    // Overwrite n_e, n_T, n_ion with functions that vary linearly along one axis

    // Density and Temperature ranges (log vals)
    const BoutReal logn_min = std::log(1e14), logn_max = std::log(1e22);
    const BoutReal logT_min = std::log(0.1), logT_max = std::log(2e4);
    const BoutReal logv_min = std::log(1), logv_max = std::log(100);

    // Linear functions for various fields that are inputs to the reaction transforms
    state["species"][heavy_reactant]["density"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logn_min, logn_max, linfunc_axis::x), &state, mesh);
    state["species"]["e"]["density"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logn_min, logn_max, linfunc_axis::y), &state, mesh);
    state["species"]["e"]["temperature"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logT_min, logT_max, linfunc_axis::z), &state, mesh);
    state["species"][atom]["velocity"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logv_min, logv_max, linfunc_axis::x), &state, mesh);
    state["species"][heavy_reactant]["temperature"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logT_min, logT_max, linfunc_axis::y), &state, mesh);
    state["species"][ion]["velocity"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logv_min, logv_max, linfunc_axis::z), &state, mesh);
    return state;
  }
};

class AmjuelHRecTest : public AmjuelReactionTest<AmjuelHydRecombinationIsotope<'h'>> {
public:
  AmjuelHRecTest()
      : AmjuelReactionTest<AmjuelHydRecombinationIsotope<'h'>>("Hrec", "h+ + e -> h",
                                                               "h+") {}
};

class AmjuelDRecTest : public AmjuelReactionTest<AmjuelHydRecombinationIsotope<'d'>> {
public:
  AmjuelDRecTest()
      : AmjuelReactionTest<AmjuelHydRecombinationIsotope<'d'>>("Drec", "d+ + e -> d",
                                                               "d+") {}
};

class AmjuelTRecTest : public AmjuelReactionTest<AmjuelHydRecombinationIsotope<'t'>> {
public:
  AmjuelTRecTest()
      : AmjuelReactionTest<AmjuelHydRecombinationIsotope<'t'>>("Trec", "t+ + e -> t",
                                                               "t+") {}
};

class AmjuelHIznTest : public AmjuelReactionTest<AmjuelHydIonisationIsotope<'h'>> {
public:
  AmjuelHIznTest()
      : AmjuelReactionTest<AmjuelHydIonisationIsotope<'h'>>("Hizn", "h + e -> h+ + 2e",
                                                            "h") {}
};

class AmjuelDIznTest : public AmjuelReactionTest<AmjuelHydIonisationIsotope<'d'>> {
public:
  AmjuelDIznTest()
      : AmjuelReactionTest<AmjuelHydIonisationIsotope<'d'>>("Dizn", "d + e -> d+ + 2e",
                                                            "d") {}
};

class AmjuelTIznTest : public AmjuelReactionTest<AmjuelHydIonisationIsotope<'t'>> {
public:
  AmjuelTIznTest()
      : AmjuelReactionTest<AmjuelHydIonisationIsotope<'t'>>("Tizn", "t + e -> t+ + 2e",
                                                            "t") {}
};

class AmjuelHeIzn01Test : public AmjuelReactionTest<AmjuelHeIonisation01> {
public:
  AmjuelHeIzn01Test()
      : AmjuelReactionTest<AmjuelHeIonisation01>("HeIzn01", "he + e -> he+ + 2e", "he") {}
};

class AmjuelHeRec10Test : public AmjuelReactionTest<AmjuelHeRecombination10> {
public:
  AmjuelHeRec10Test()
      : AmjuelReactionTest<AmjuelHeRecombination10>("HeRec10", "he+ + e -> he", "he+") {}
};

#endif