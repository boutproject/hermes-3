# Relaxing radial potential boundaries

This nonperiodic X-Z slab evolves two ions with unequal mass-to-charge ratios
and a common temperature. Density/composition and electron pressure vary in
X-Z. Ion parallel momenta are nonzero, so polarisation transports both enthalpy
and parallel kinetic energy. Y is periodic and the initial fields are uniform
in Y; Z is periodic. Two Y guards support the parallel transport stencils.

The potential boundary midpoint relaxes towards the adjacent interior Z
average on a 1 microsecond timescale. Its state is carried by the physical
boundary cells of solver `Vort`, while ordinary output `Vort` and stencil
inputs contain physical vorticity. The auxiliary current potential has
homogeneous midpoint Dirichlet DC and radial Neumann AC boundaries.
The full current operator permits DC and mixed-metric normal current, but
particle, momentum, enthalpy, parallel kinetic and
polarisation electrical fluxes are closed at physical radial faces.

The independent Python audit reconstructs inventories, the full current solve,
each closed polarisation contribution, midpoint relaxation rates and all
supported source/energy face fluxes. It includes any Reynolds/electrical
boundary power rather than assuming that all boundary powers vanish. External
sources are prescribed zero and validated against the input.

Run with the workspace environment:

```sh
.venv/bin/python tests/integrated/2D-ion-energy-relax-boundary/runtest --exe build/hermes-3
```

The CI profile checks three output cadences, a diagnostics-off pair, an
optional-diagnostics-stripped audit, two-rank radial and poloidal decompositions and a
checkpoint/restart against an uninterrupted run. MPI uses `MPIEXEC` when set,
otherwise `mpiexec -n`; the executable and results directory are configurable.
Each run writes isolated outputs, logs and JSON assertions. `--profile smoke`
uses the two finest cadences while retaining boundary and restart assertions.

Diagnostic comparisons use tighter tolerances than the restart comparison.
The auxiliary chi audit reproduces the fixed spectral corrections and reports
its FV residual separately from reconstruction agreement. This smooth fixture
uses a relative residual budget of `1e-3`; the exact FV reference must still
satisfy the original forward tolerance. Boundary closure and power balance
retain their independent strict checks.
A restart resets CVODE history, so its comparison permits differences within
the configured time solver accuracy (`rtol=1e-9`), using `1e-7` relative and
`1e-8` absolute history tolerances. The midpoint ODE audit independently
checks relaxation in every uninterrupted, MPI and resumed run.
Additional runs tighten CVODE tolerances by ten and halve the relaxation
timescale. Tightened solver histories must agree within a `1e-7` normalised
absolute-plus-relative accuracy bound, including small vorticity and
potential values where local relative error alone is misleading. Reports
retain the maximum error for each field. Halving the timescale must double
the initial reconstructed rate.
Both runs retain the full face and power audits.

The test uses the same physical vorticity boundary configuration (`neumann`)
for stencil reconstruction, independent of the solver storage's no-op
operators. Replica consistency, core averages, mixed-metric face closure and
solver boundary operations also have C++ unit coverage.
