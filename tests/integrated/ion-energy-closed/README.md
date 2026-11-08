# Closed two-ion energy balance

A periodic 16 x 4 x 16 slab has deuterium, a doubly charged impurity,
quasineutral zero-current electrons, curvature and viscosity. Composition and
parallel momenta vary in X/Z, and the momenta also vary in Y to drive finite
parallel-energy exchange. Initial total ion energy includes the prescribed
parallel kinetic energy so the thermal reservoir stays positive. Diagnostics on/off runs check energy, species
particles, charge, recovery and reconstructed diagnostics at a fixed `1e-10`
floor. Curvature exchange and both kinetic-energy changes must exceed `1e-4`
of the initial energy, with tenfold margins in the calibrated profile.

```sh
python tests/integrated/ion-energy-closed/runtest --exe build/hermes-3
```

`criteria.json` records limits and observations. Isolated `result.json` files
retain effective limits, summaries, input hashes and execution times. There
are no sources, sheaths or parallel current. Physics failures remain enabled.

Known failure: both diagnostics settings show relative energy drift
`1.426e-6` against `1e-10` (cause unknown). The check remains enabled.
