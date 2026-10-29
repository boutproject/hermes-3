# Interchange publication study

`publication_scan.py` runs matched pressure and total-ion-energy formulations
and retains the inputs, time-zero restart files, raw dumps, logs and scalar
analysis. The original `data/BOUT.inp` and `runtest` remain the historical
pressure-sign regression.

Use the repository environment from the repository root:

```sh
.venv/bin/python tests/integrated/2D-energy/publication_scan.py \
  --exe build/hermes-3 --phi-boundary both --accuracy-pilot \
  --output-root /path/to/results

# Bounded validation: 8² and 16², fixed and relaxing boundaries, viscosity=0.
.venv/bin/python tests/integrated/2D-energy/publication_scan.py --smoke

# Recreate plots and tables without executing Hermes.
.venv/bin/python tests/integrated/2D-energy/publication_scan.py \
  --analyse-only /path/to/results/interchange-TIMESTAMP --latex
```

The default is four physical grids (32², 64², 128², 256²), fixed phi boundaries,
viscosity 1 m²/s, duration 20000 normalised time units and 200 output intervals.
`--phi-boundary relax` selects solver-managed relaxation;
`--phi-boundary both` runs both policies. The relaxation timescale defaults to
1e-4 seconds and is independent of resolution. `--levels`, `--base-nx`,
`--base-nz`, `--refinement`, `--nproc`, `--end-time`, `--outputs`, `--exe` and
`--timeout` configure the simulation matrix. With singleton Y, MPI decomposes
only X; each rank needs at least four physical X cells on the coarsest grid.
The launcher uses the installed `mpiexec` and inherits its environment.

Four snapshot times default to fractions 0, 1/4, 1/2 and 1 of the duration;
the convergence time defaults to half the duration. Explicit
`--snapshot-times T0 T1 T2 T3` and `--convergence-time T` must be represented
exactly on the output axis. The runner refines the requested regular cadence
when necessary, then checks all actual saved times, including the first record
at zero. `--field` selects `ion_pressure`, `density` or `vorticity`.

Every run gets its own generated `BOUT.inp` from `publication-input.inp`.
The study manifest records configuration, input/executable/runner hashes,
hashes of separately built Hermes/BOUT libraries, checkout revisions and dirty
state, Python environment and MPI layout.
`--resume STUDY` requires matching simulation arguments, executable and source
hashes; repeat the original flags (including `--smoke`, if used). The level
count can stay the same or increase; increasing it appends finer grids.
After updating the runner, add `--allow-runner-change` to accept and record its
new hash; configuration, input/template and executable/library checks still
apply, and the original runner hash is retained in the manifest.
Completed runs are reused; incomplete outputs are archived before a fresh attempt. Runtime limits and
interrupts stop the local launch group and retain an explicit failure status.
`--analyse-only` uses the retained configuration and rejects changed inputs,
missing records and irregular output axes.

A nonzero process exit can occur after the final output, for example during
MPI finalisation. The runner checks every rank's complete expected time axis
and finite required physical fields in the final record. If these checks pass,
it accepts the outputs and records `outputs_complete=true` in `status.json`,
preserving the nonzero exit code and printing a warning. Startup recovery also
requires readable restart files for every rank. Resume and reanalysis perform
the same check for previously failed runs. Incomplete or unreadable outputs
remain failures; a keyboard interrupt stops the current invocation even when
its outputs are complete.

## Extending and distributing a study

To add resolution levels without repeating all the original simulation flags,
use `--extend`. The number of levels is the **new total**, not the number to
add. All saved simulation settings and the original executable path are
inherited; only an increased level count is accepted.

```sh
.venv/bin/python tests/integrated/2D-energy/publication_scan.py \
  --extend /path/to/results/interchange-TIMESTAMP --levels 5 \
  --allow-runner-change
```

The manifest retains existing cases and records the appended grids. Completed
cases are reused; remaining cases run normally. Once the whole matrix is
complete, figures and tables are regenerated, using each formulation's new
finest grid as its convergence reference. Decreasing the level count or
changing physical settings is rejected. The runner-version override is needed
only when the script has changed since the study was created or last resumed.

For remote execution, prepare matched initial conditions locally:

```sh
.venv/bin/python tests/integrated/2D-energy/publication_scan.py \
  --extend /path/to/results/interchange-TIMESTAMP --levels 5 \
  --allow-runner-change --prepare-only
```

`--prepare-only` also works with a new study or `--resume`. It runs the
zero-output startup simulations needed to generate native restart metadata,
then writes matched time-zero restarts and a standalone `run-remote.py` into
each pending case. It stops before the full simulations and analysis.
`--cases CASE [CASE ...]` selects particular directory names, such as
`fixed-pressure-128x128` and `fixed-total_energy-128x128`, for preparation or
execution. Completed cases are skipped. Prepared restarts are retained and
can also be used by an ordinary local resume.

Copy these files from each prepared case to a directory on the other machine:

- `BOUT.inp`;
- all `BOUT.restart.<rank>.nc` files;
- `initial-state.nc`, `remote-case.json` and `run-remote.py`.

The nested `startup/` directory is not needed remotely. Keep the original
`initial-state.nc` and `remote-case.json` locally. The remote simulation needs
the same Hermes/BOUT source revisions and MPI rank count, but can use a build
for a different architecture. The helper needs only Python's standard library:

```sh
python3 /remote/cases/fixed-pressure-128x128/run-remote.py \
  --exe /remote/hermes/build/hermes-3 --source-root /remote/hermes
```

The helper verifies the prepared file hashes and uses the saved MPI process
count. `--mpiexec PATH` selects a launcher; `--timeout SECONDS` sets an optional
wall-clock limit. MPI inherits the remote environment. The helper retains
stdout/stderr in `run.log` and writes `remote-execution.json` with the exit
code, hostname, executable/library hashes and checkout revisions. It refuses
to overwrite an existing attempt. Copy `BOUT.inp`, every rank's
`BOUT.dmp.<rank>.nc`, `BOUT.settings`, logs, and `remote-execution.json` back under their original
case directory names in a return directory, then import them:

```sh
.venv/bin/python tests/integrated/2D-energy/publication_scan.py \
  --import-runs /path/to/results/interchange-TIMESTAMP \
  --import-from /path/to/returned-cases \
  --cases fixed-pressure-128x128 fixed-total_energy-128x128
```

Alternatively copy outputs directly into the original local case directories
and omit `--import-from`. No local Hermes executable or MPI launcher is needed
for importing. Without `--cases`, import scans pending cases with returned
dumps. Existing completed cases are reused. Each imported case must have the
original input, every MPI rank's full time axis, finite final state, the
expected grid and physical box, and the matched initial pressure, density,
potential, vorticity and (for total energy) `Eion`. External outputs are staged
and checked before replacing a local failed attempt or writing completion
status. Known source-revision mismatches or a build changed during execution
are rejected. A nonzero remote exit is retained and can be accepted when the
outputs pass these checks.

Manual remote launches are supported too: return the input, dumps and logs
and import them the same way. Return `BOUT.settings` when available: its
resolved values must agree with every supplied input option. If it is absent,
loading falls back to `BOUT.inp`, and the import records that only supplied
options were available. If `remote-execution.json` is absent, the script
records the exit code and remote build provenance as unknown and prints a
warning. It does not invent a successful exit code. Use the helper when
collecting publication evidence.

Import registers completion without running remaining simulations or doing
analysis. Use `--resume` to complete the matrix, or `--analyse-only` when all
cases have finished. Copied dumps alone are not treated as completed by a
resume; import them first.

## Physical setup and startup

Both models evolve hydrogen density, impose quasineutral electrons, and evolve
electron pressure. The pressure formulation also evolves ion pressure and uses
full Reynolds stress and legacy polarisation advection. The total-energy model
evolves `Eion` and recovers the common ion temperature. There is no parallel
flow, conduction, source or hyperdiffusion. The box is 5 cm × 5 cm, with
R=1.5 m, B=0.35 T, Nnorm=2e18 m⁻³ and Tnorm=5 eV. The template supplies consistent SI metric pairs explicitly and uses
`normalise_metric=true`; the legacy `recalculate_metric` path mixes B units
when constructing its Z metric for this setup. Guard widths stay fixed when
physical grid counts are refined.

Initial density is uniform and both species have
`P=1+0.1*sin(2*pi*x-z)`. The shared physical potential is sampled from the
analytic continuum ion-only zero-vorticity solution: `Phi` is the harmonic
extension of the pressure boundary values and `phi=Phi-P`. Apply the same
homogeneous Neumann pressure guards and zero phi midpoints in both models.
Each model's initial `Vort` is obtained by applying its own forward operator
to the sampled generalised potential. For total energy, initialise
`Eion=1.5*Pi+Kperp_FV`; parallel kinetic energy is zero. A zero-output seed run
provides native MPI restart metadata; generated restarts replace its fields
and reset time and iteration to zero. `initial-state.nc` retains physical fields
separately from the restart boundary-state encoding. The first saved record is
checked against those arrays.

This is a comparison of two specified formulations. The pressure model uses
an approximate potential-flow polarisation transport and includes electron
pressure in its mass/charge-weighted generalised potential. The total-energy
model uses an ion-current partition and ion-only pressure potential. Their
native inversion and viscous heating discretisations also differ. The pressure
sine is incompatible with homogeneous radial Neumann derivatives at startup;
boundary-adjacent physical cells are retained in all primary error norms.
These differences are recorded in every manifest.

## Energy and convergence analysis

The common physical estimator uses ion-only `psi`, the same finite-volume
`calculateKperp` reconstruction in both models, ion/electron thermal energy and
zero parallel kinetic energy. It is distinct from each model's native energy.
For total energy, native plasma energy is saved `Eion` plus electron thermal
energy. For pressure, the native quadratic energy uses FD X and spectral Z,
including the Nyquist mode and radial half-face terms. All energies are in
joules for the represented box, integrated over physical cells only. No
additional toroidal multiplier is applied.

Total-energy face powers are reconstructed independently, including Reynolds,
electrical and viscous transport. Polarisation transport is closed at radial
faces while polarisation current remains in the vorticity equation. The native
residual is `E(t)-E(0)+integral(Pout dt)` because prescribed sources are zero.
No extra pressure-model boundary-relaxation estimate is added to this balance.
For production inventories, closed normal polarisation/ExB faces and periodic Z
allow global outward power to be evaluated without an auxiliary solve. This
shortcut has an explicit singleton-Y, zero-parallel-flow, closed-species-face
contract; it does not replace locally nonzero transport divergences with zero.
Diagnostic runs reconstruct the full operators to check the independent audit.
Pressure runs retain raw native/common drift and the independently reconstructed
viscous heating defect. Their corrected residual is **unavailable**: full legacy
polarisation and time-dependent `Phi` boundary work are not covered. An absent
power is never substituted with zero.

`--diagnose` enables optional cross-checks; production defaults to false.
Inventories and powers use ordinary state fields and metrics. Analysis streams
time records and retains scalar histories, loading state arrays for the four
images and declared convergence time. Saved phi guards also provide boundary
midpoint histories and a trapezoid residual of the relaxation ODE.

Global convergence uses common `Kperp` at the declared time, with signed adjacent
differences, finest-reference differences and observed orders. Orders are
omitted when adjacent differences change sign or reach roundoff. Each
formulation uses its own finest reference. Pressure/density use conservative
volume-overlap restriction; periodic Z wrapping handles the shifted Fourier
cell faces. Vorticity uses Fourier Z and cubic X point sampling. Local errors
are volume-weighted RMS (`L2`) and maximum (`Linf`), plus relative norms with
explicit zero-reference handling. The finest-reference error is unavailable,
rather than a zero on logarithmic axes.

`--accuracy-pilot` repeats the two finest grids with tighter time/potential
solver tolerances and half the maximum timestep, then with twice the output
cadence. `accuracy-pilot.json` records observable and field changes relative to
spatial differences, with a target below 10%, and compares cumulative energy
residuals. Freeze the final settings after inspecting this evidence. Failure to
self-converge or a remaining energy drift is a reported numerical result.
A short smoke run does not validate temporal error for the full-duration study.

## Products

Each boundary policy produces four PDFs:

- `state-*.pdf`: both formulations at four common times on their finest grid,
  with shared colour limits, physical coordinates and rasterised image artists.
- `energy-*.pdf`: thermal and perpendicular kinetic changes, common and native
  total drift, available corrected residuals and native/common discrepancy.
  Panel labels call these "shared FV energy" and "model energy"; the finest
  grid uses solid traces and coarser grids use distinct dash patterns.
- `global-convergence-*.pdf`: common `Kperp`, adjacent differences and valid
  observed orders at the specified time.
- `local-error-*.pdf`: absolute RMS and maximum error against each formulation's
  own finest reference.

`energy-table` and `convergence-table` are printed and saved as CSV and JSON.
`--latex` adds tables from the same numeric data. Field units are retained in
the convergence data. Each run also retains `inventory.nc` and `analysis.json`,
including diagnostic comparisons and boundary/viscous audit information.
`figure-config.json` records plotting settings; `--style FILE` accepts a JSON
mapping of Matplotlib rcParams. Reanalysis reuses the saved settings; an explicit
`--style` overrides them. Fonts are embedded in the PDF and traces/labels remain
vector graphics.
