#include "test_adas_reactions.hxx"
#include "test_amjuel_reactions.hxx"

// H isotopes ionization
TEST_F(AmjuelHIznTest, SourcesRegression) { sources_regression_test(); }
TEST_F(AmjuelDIznTest, SourcesRegression) { sources_regression_test(); }
TEST_F(AmjuelTIznTest, SourcesRegression) { sources_regression_test(); }

// H isotopes recombination
TEST_F(AmjuelHRecTest, SourcesRegression) { sources_regression_test(); }
TEST_F(AmjuelDRecTest, SourcesRegression) { sources_regression_test(); }
TEST_F(AmjuelTRecTest, SourcesRegression) { sources_regression_test(); }

// H isotopes CX (non-exhaustive)
TEST_F(HHpCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(DDpCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(TTpCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(HDpCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(THpCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(DTpCXTest, SourcesRegression) { sources_regression_test(); }

// He ionization
TEST_F(AmjuelHeIzn01Test, SourcesRegression) { sources_regression_test(); }

// He recombination
TEST_F(AmjuelHeRec10Test, SourcesRegression) { sources_regression_test(); }

// C ions ionization (non-exhaustive)
TEST_F(ADASCIznTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASCpIznTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASC5pIznTest, SourcesRegression) { sources_regression_test(); }

// C ions recombination (non-exhaustive)
TEST_F(ADASCRecTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASCpRecTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASC5pRecTest, SourcesRegression) { sources_regression_test(); }

// C ions CX  (non-exhaustive)
TEST_F(ADASCpHCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASCp3DCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASCp5TCXTest, SourcesRegression) { sources_regression_test(); }

// Li ions ionization
TEST_F(ADASLiIznTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASLipIznTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASLi2pIznTest, SourcesRegression) { sources_regression_test(); }

// Li ions recombination
TEST_F(ADASLiRecTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASLipRecTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASLi2pRecTest, SourcesRegression) { sources_regression_test(); }

// Li ions CX (non-exhaustive)
TEST_F(ADASLipHCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASLip2DCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASLip3TCXTest, SourcesRegression) { sources_regression_test(); }

// Ne ions ionization (non-exhaustive)
TEST_F(ADASNeIznTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASNepIznTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASNe9pIznTest, SourcesRegression) { sources_regression_test(); }

// Ne ions recombination (non-exhaustive)
TEST_F(ADASNeRecTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASNepRecTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASNe9pRecTest, SourcesRegression) { sources_regression_test(); }

// Ne ions CX  (non-exhaustive)
TEST_F(ADASNepHCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASNep5DCXTest, SourcesRegression) { sources_regression_test(); }
TEST_F(ADASNep9TCXTest, SourcesRegression) { sources_regression_test(); }