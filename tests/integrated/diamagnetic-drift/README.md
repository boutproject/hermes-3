# Diamagnetic drift on a single-null tokamak

The COMPASS 36x48 grid supplies the metric and curvature. A small model
evolves only production diamagnetic density, pressure and momentum sources,
with unequal target temperature guards. It checks volume-weighted conservation
and zero source integrals to relative 1e-12, nonzero evolution, closed radial
flux deposition and uniform core-ring fields/sources. Core averaging preserves
the ring integral and leaves other cells unchanged relative to the production
`closed` policy. No discrete transport operator is copied into the test.

Native fields and source diagnostics must agree to relative 1e-10 on layouts
(NXPE,NYPE)=(1,6), (1,12), (2,6), for MZ=1 and MZ=4. The grid's eight-cell
divertor legs require at least six Y processors. Unit properties separately
cover gradient/bracket forms, charge signs, skew metrics and sheath closure.

Build `diamagnetic_drift_smoke`, then run `ctest -R '^diamagnetic-drift$'`
with boutdata on PATH. CTest supplies executable/grid paths and a timeout.
The runner uses `launch_safe` and accepts `--exe`, `--grid` and
`--output-root` (a fresh directory). Inputs, native outputs and logs are kept.
