# Physical sources and ion energy conversion

A uniform periodic two-ion slab combines electron pressure heating, positive
and negative density sources on moving ions, and Braginskii collisional heat
exchange. The shared analysis reconstructs source powers from resolved options
and saved state, checks energy and particle balances, and compares two output
cadences to estimate trapezoidal sampling error. Diagnostics on/off runs must
produce the same final state.

```sh
python tests/integrated/ion-energy-sources/runtest --exe build/hermes-3
```

`result.json` contains source terms, coverage, cadence estimates, hashes, logs
and all checks. `criteria.json` records limits and observations.
No sheaths or parallel current are present. Friction is omitted because the
zero-current electron closure has no evolved electron kinetic reservoir for
its mechanical work. Existing C++ injection/removal unit tests cover momentum
source conversion. Collisional heat exchange is an internal transfer.
