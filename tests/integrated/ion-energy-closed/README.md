# Closed two-ion energy balance

A periodic 16 x 4 x 16 slab has deuterium, a doubly charged impurity,
quasineutral zero-current electrons, constant curvature and viscosity. Composition,
parallel momenta and vorticity vary across the drift plane. Two identical runs
with diagnostics enabled and disabled check energy, species particles, charge,
recovery consistency and nontrivial changes in the kinetic and thermal reservoirs.

```sh
python tests/integrated/ion-energy-closed/runtest --exe build/hermes-3
```

The shared analysis writes inputs, hashes, logs, summaries and check values to
an isolated `result.json`. `criteria.json` holds physical limits and measured
values. There are no sources, sheaths or parallel current; the detailed Reynolds
flux is not an acceptance condition.
