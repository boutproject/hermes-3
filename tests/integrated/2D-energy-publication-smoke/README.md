# Publication runner smoke test

Run the [interchange publication study](../2D-energy/publication-study.md) on
8² and 16² physical grids, in both pressure and total-ion-energy formulations,
with fixed and relaxing phi boundaries. Duration is 12 normalised time units,
four output intervals, and viscosity is explicitly **zero**. This is a bounded
pipeline/startup regression, not evidence of full-duration spatial convergence.

```sh
.venv/bin/python tests/integrated/2D-energy-publication-smoke/runtest \
  --exe build/hermes-3
```

Inputs are generated independently for each run from the shared publication
input template and retained with its manifest. This test checks matching physical
time-zero fields, exact output times, fixed potential midpoints, positive and
resolution-independent represented volumes, reconstructed total-energy
partitions and a relative native total-energy drift below 1e-7. It verifies all
eight PDF products, CSV/JSON tables, and reanalysis from retained outputs.
Pressure corrected residuals must be marked unavailable until legacy boundary
and polarisation power coverage is complete. Results remain in the printed
output directory; `--output-root` chooses a persistent parent directory.
