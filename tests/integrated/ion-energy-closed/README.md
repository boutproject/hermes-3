# Closed two-ion energy balance

A periodic 16 x 4 x 16 slab has deuterium, a doubly charged impurity,
quasineutral zero-current electrons, curvature and viscosity. Composition and
parallel momenta vary in X/Z. All fields are invariant in Y. Initial total ion
energy includes the prescribed parallel kinetic energy so the thermal reservoir
stays positive. Diagnostics on/off runs check energy, species
particles, charge, recovery and reconstructed diagnostics at a fixed `1e-10`
floor. Curvature exchange and perpendicular kinetic-energy change must exceed
`1e-4` of initial energy; parallel kinetic-energy change must exceed `5e-5`.
The parallel minimum was revised with maintainer approval for the Y-invariant
benchmark.

```sh
python tests/integrated/ion-energy-closed/runtest --exe build/hermes-3
```

`criteria.json` records limits and observations. Isolated `result.json` files
retain effective limits, summaries, input hashes and execution times. There
are no sources, sheaths or parallel current. Physics failures remain enabled.

The previous input also modulated momenta in Y. It showed relative energy drift
`1.42336e-6` against `1e-10`. This is an electron-closure limitation:
`zero_current` sets velocity from ion current without enforcing parallel force
balance or evolving electron kinetic energy. The unpaired work is
`Ve Grad_par(Pe) + Je Grad_par(phi)`, which is not forced to vanish by zero total
current. The strict conservation benchmark therefore uses Y-invariant fields;
parallel kinetic energy still changes through perpendicular transport. The
original failed inputs/outputs are retained with the polarisation-drift review.
The analysis reports this work separately from heating and omits a closed
residual when it is nonzero or cannot be reconstructed. Parallel force-balance
and evolved-electron tests belong to the 3D current/MMS stage of that plan.

The reviewed X-Z input conserves energy with measured relative drift
`1.37519e-16`; its parallel kinetic energy change is `8.53316e-5`, above the
revised `5e-5` minimum. The maintainer accepted checkpoint 0 with this revision.
Historical Y-varying and reviewed X-Z measurements are distinguished in
`criteria.json`; the original failing evidence remains intact.
