#include "test_amjuel_reactions.hxx"

// H isotopes ionization
TEST_F(AmjuelHIznTest, SourcesRegression) { sources_regression_test(); }
TEST_F(AmjuelDIznTest, SourcesRegression) { sources_regression_test(); }
TEST_F(AmjuelTIznTest, SourcesRegression) { sources_regression_test(); }

// H isotopes recombination
TEST_F(AmjuelHRecTest, SourcesRegression) { sources_regression_test(); }
TEST_F(AmjuelDRecTest, SourcesRegression) { sources_regression_test(); }
TEST_F(AmjuelTRecTest, SourcesRegression) { sources_regression_test(); }

// He ionization
TEST_F(AmjuelHeIzn01Test, SourcesRegression) { sources_regression_test(); }

// He recombination
TEST_F(AmjuelHeRec10Test, SourcesRegression) { sources_regression_test(); }