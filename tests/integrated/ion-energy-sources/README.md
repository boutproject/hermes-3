# Physical sources and ion energy conversion

A periodic two-ion slab combines electron heating, density injection/removal
on moving ions and Braginskii heat exchange. Analysis reconstructs source
powers from resolved options and saved state, and uses two output cadences to
estimate trapezoidal sampling error. After subtracting that estimate, the
energy residual uses `max(1e-10, rtol_factor * solver:rtol)`. Particle balances,
recovery and diagnostics on/off state differences use a fixed `1e-10` limit.

```sh
python tests/integrated/ion-energy-sources/runtest --exe build/hermes-3
```

`result.json` contains effective limits, source coverage, cadence estimates,
input hashes and execution times; `criteria.json` records observations.
There are no sheaths or parallel current. Friction is omitted because the
zero-current electron closure has no evolved kinetic reservoir for its work.
Heat exchange is an internal transfer; C++ tests cover momentum-source conversion.
