#ifndef TEST_AMJUEL_REACTIONS_H__
#define TEST_AMJUEL_REACTIONS_H__

#include "test_reactions.hxx"

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
  AmjuelReactionTest(std::string lbl, std::string reaction_str, std::string sp_in)
      : ReactionTest<RTYPE>(lbl, reaction_str), sp_in(sp_in) {}

protected:
  std::string sp_in;

  Options generate_state() override {
    std::string atom = this->sp_in.substr(0, 1);
    std::string ion = atom + "+";

    const std::map<std::string, BoutReal> sp_masses = {
        {"h", 1.0}, {"d", 2.0}, {"t", 3.0}, {"e", 1. / 1836}};
    const BoutReal atom_charge = 0.0;
    const BoutReal ion_charge = 1.0;
    std::string comp_name("test" + this->lbl);
    Options state{{comp_name, {{"type", this->reaction_str}}},
                  {"units", {{"eV", 1.0}, {"inv_meters_cubed", 1.0}, {"seconds", 1.0}}},
                  {"species",
                   {{"e",
                     {{"AA", sp_masses.at("e")},
                      {"density", 1.0},
                      {"temperature", 1.0},
                      {"velocity", 1.0}}},
                    {atom,
                     {{"AA", sp_masses.at(atom)},
                      {"charge", atom_charge},
                      {"density", 1.0},
                      {"temperature", 1.0},
                      {"velocity", 1.0}}},
                    {ion,
                     {{"AA", sp_masses.at(atom)},
                      {"charge", ion_charge},
                      {"density", 1.0},
                      {"temperature", 1.0},
                      {"velocity", 1.0}}}}}};

    // Overwrite n_e, n_T, n_ion with functions that vary linearly along one axis

    // Density and Temperature ranges (log vals)
    const BoutReal logn_min = std::log(1e14), logn_max = std::log(1e22);
    const BoutReal logT_min = std::log(0.1), logT_max = std::log(2e4);

    constexpr BoutReal xmin = 0, xmax = 2;
    constexpr BoutReal ymin = 0, ymax = 25.1327412287;
    constexpr BoutReal zmin = 0, zmax = 5.38558740615;

    // Linear functions for input atom/ion density, e density, e temperature
    state["species"][sp_in]["density"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logn_min, logn_max, "x", xmin, xmax), &state, mesh);
    state["species"]["e"]["density"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logn_min, logn_max, "y", ymin, ymax), &state, mesh);
    state["species"]["e"]["temperature"] = FieldFactory::get()->create3D(
        this->gen_lin_field_str(logT_min, logT_max, "z", zmin, zmax), &state, mesh);

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

#endif