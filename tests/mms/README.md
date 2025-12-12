MMS tests of differential operators used in Hermes-3
====================================================

These tests are designed to test that for any given
differential operator of two arguments `L(a,f)` 
returns the expected result for a known `a`, `f`, and
contravariant metric coefficients, i.e., we compute
a symbolic result

    S = L(a,f)

from prescribed symbolic expressions using `sympy`,
and we check that the numerical result for `S` is equal to the
symbolic one up to the numerical error.

We carry out two tests, one on operators which are designed
for "orthogonal" metrics (`orthogonal_test.py`), and another for operators which
are generalised to fully nonorthogonal metrics (`nonorthogonal_test.py`).

For each operator, we plot the numerical error, showing
how it declines compared to the expected convergence order
of 2 for the conservative finite difference methods.
These plots can be made interactively, or are otherwise plotted
to `fig_{i}.png` for the `ith` operator.

The tests can be extended to include new operators by including
new operators in `const auto differential_operators` in `main.cxx`
and in `"differential_operator_list":` from the `test_input` dictionary
in the python test scripts `orthogonal_test.py` or `nonorthogonal_test.py`.
