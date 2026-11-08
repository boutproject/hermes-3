# Interchange publication study

`publication_scan.py` runs matched pressure and total-ion-energy formulations
and retains the inputs, time-zero restart files, raw dumps, logs and scalar
analysis. Energy drift and spatial self-convergence are separate measurements:
conserving total energy does not establish spatial accuracy. The original
`data/BOUT.inp` and `runtest` remain the historical pressure-sign regression.

Use the repository `.venv` with editable xhermes and the numerical/plotting
dependencies. From the repository root, install these into an existing
environment:

```sh
uv pip install --python python -e external/xhermes \
  -r tools/energy_analysis/requirements.txt
```

Run the study from the repository root:

```sh
python examples/tokamak-2D-driftplane/2D-energy/publication_scan.py \
  --exe build/hermes-3 --phi-boundary both --accuracy-pilot \
  --output-root /path/to/results

# Bounded validation: 8² and 16², fixed and relaxing boundaries, viscosity=0.
python examples/tokamak-2D-driftplane/2D-energy/publication_scan.py --smoke

# Recreate plots and tables without executing Hermes.
python examples/tokamak-2D-driftplane/2D-energy/publication_scan.py \
  --analyse-only /path/to/results/interchange-TIMESTAMP --latex

# Export each energy-history panel as a standalone PDF.
python examples/tokamak-2D-driftplane/2D-energy/publication_scan.py \
  --analyse-only /path/to/results/interchange-TIMESTAMP --energy-format pdf-panels
```

## Run configuration

The default matrix has eight runs: both formulations on four physical grids
(32², 64², 128², 256²), with fixed phi boundaries. `--phi-boundary both`
produces sixteen runs. Cases run sequentially and each owns its inputs and
outputs. The publication matrix is an explicit longer study; `--smoke` is
the bounded, inviscid validation profile.

| Option | Default | Meaning |
| --- | --- | --- |
| `--levels` | `4` | Number of resolutions; at least two |
| `--base-nx`, `--base-nz` | `32`, `32` | Physical X and Z cells on the coarsest grid; at least four |
| `--refinement` | `2` | Integer refinement ratio in both directions; at least two |
| `--phi-boundary` | `fixed` | `fixed`, `relax`, or `both` |
| `--phi-boundary-timescale` | `1e-4` | Pressure-model wall relaxation time [s] |
| `--phi-boundary-velocity` | `1000` | Total-energy wall relaxation velocity [m/s] |
| `--viscosity` | `1` | Common physical kinematic viscosity [m²/s] |
| `--end-time` | `20000` | Duration in normalised time, `t*omega_ci` |
| `--outputs` | `200` | Regular output intervals, plus the initial record; at least four |
| `--snapshot-times` | `0, T/4, T/2, T` | Four normalised snapshot times, where `T` is the duration |
| `--convergence-time` | `T/2` | Normalised time for global and local convergence |
| `--field` | `ion_pressure` | Local convergence field; also `density` or `vorticity` |
| `--nproc` | `1` | MPI ranks, decomposing X only |
| `--timeout` | `7200` | Wall-clock limit [s] for each startup or simulation launch |
| `--exe` | `hermes-3` or `build/hermes-3` | Repository executable, preferring the former when present |
| `--output-root` | `publication-results` beside the script | Parent of a new `interchange-TIMESTAMP` study |

`--phi-boundary relax` selects solver-managed relaxation in both formulations.
The pressure model uses `(mean(phi_adjacent)-phi_b)/tau`, with
`tau=--phi-boundary-timescale`. Total ion energy selects
`radial_boundary=no_slip` and uses `dphi_b/dt=kappa*dphi/dn`, where `n` points
into the domain and `kappa=--phi-boundary-velocity`. Both start with zero
potential midpoints. Fixed cases hold those midpoints constant. Each boundary
parameter stays fixed across resolutions; the two laws do not imply equal
relaxation strength.

With singleton Y, MPI decomposes only X. The coarsest physical X count must
be divisible by `--nproc`, with at least four physical X cells per rank.
The launcher uses the installed `mpiexec` and inherits its environment.

Time integration and total-energy potential recovery have separate controls:

| Option | Default | Meaning |
| --- | --- | --- |
| `--atol` | `1e-10` | CVODE absolute tolerance |
| `--rtol` | `1e-5` | CVODE relative tolerance |
| `--max-timestep` | `50` | Maximum internal timestep in normalised units |
| `--phi-atol` | `1e-3` | Total-energy potential iteration scale [eV] |
| `--phi-rtol` | `1e-5` | Total-energy potential recovery relative tolerance |
| `--vort-atol` | `1e-10` | Total-energy normalised vorticity residual tolerance |

`--phi-atol` is the scale in the relative potential-update denominator;
it is not a standalone bound on the recovered field error. Use the accuracy
pilot described below to compare solver and output errors with spatial errors.

Total-ion-energy inputs select the wall model through `radial_boundary`;
Laplacian boundary flags are supplied only to the pressure formulation. When
running the script from `build/tests/integrated/2D-energy`, refresh the copied
script and template after source changes with `cmake -S . -B build` from the
repository root.

Four snapshot times default to fractions 0, 1/4, 1/2 and 1 of the duration;
the convergence time defaults to half the duration. Explicit
`--snapshot-times T0 T1 T2 T3` and `--convergence-time T` must be represented
exactly on the output axis. The runner refines the requested regular cadence
when necessary, then checks all actual saved times, including the first record
at zero. For example, `--outputs 50` becomes 100 intervals when the default
quarter-duration snapshot is required. The refined cadence is recorded in
the manifest. `--field` selects `ion_pressure`, `density` or `vorticity` for the
convergence analysis. With the default `ion_pressure`, the upper snapshot row
shows ion internal energy density `3*Pi/2` from pressure evolution; other
choices plot the selected field.
The lower row shows perpendicular kinetic energy density from total ion energy.

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
python examples/tokamak-2D-driftplane/2D-energy/publication_scan.py \
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
python examples/tokamak-2D-driftplane/2D-energy/publication_scan.py \
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
python examples/tokamak-2D-driftplane/2D-energy/publication_scan.py \
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
flow, conduction, source or hyperdiffusion. Both use reference atomic mass 1.
The box is 5 cm × 5 cm, with R=1.5 m, B=0.35 T, Nnorm=2e18 m⁻³ and Tnorm=5 eV.
The template supplies consistent SI metric pairs explicitly and uses
`normalise_metric=true`; the legacy `recalculate_metric` path mixes B units
when constructing its Z metric for this setup. Guard widths stay fixed when
physical grid counts are refined: `MXG=2`, `MYG=0`, `mesh:ny=1` and
`mesh:nx=nx_physical+4`. Thus the default X mesh sizes are 36, 68, 132 and 260.
X has physical walls, Z is periodic, and the constant curvature is
`bxcvz=1/Rxy^2`. Cell widths are recomputed at each resolution so the box and
perturbation wavelengths remain fixed.

Initial density is uniform and both species have
`P=1+0.1*sin(2*pi*x)*sin(2*pi*x-z)`, with `x` between 0 and 1. The perturbation
vanishes at both X faces, making the initial pressure constant in Z there.
The shared physical potential is sampled from the analytic continuum ion-only
zero-vorticity solution: the harmonic extension is `Phi=1`, so `phi=1-P`.
Both models have zero phi midpoints. The pressure formulation uses homogeneous Neumann
ion-pressure guards; total ion energy imposes mean-Neumann, zero-AC-midpoint
guards on ion pressure and its pressure potential. Electron pressure retains
its species Neumann guards. These guard conditions leave the initial physical
density, pressures and potential matched.
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
native inversion and viscous heating discretisations also differ. Total energy
uses conservative face-based viscous work and heating; the pressure model's
cell-centred heating defect is reconstructed and reported. The tapered
pressure perturbation has constant X-face values, but its radial derivative
need not satisfy the pressure formulation's homogeneous Neumann condition.
Boundary-adjacent physical cells are retained in all primary error norms.
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
differences, finest-reference differences and observed orders. With physical
spacing `h`, refinement ratio `r`, and `Q_h=Kperp_common(t_convergence)`,

$$
\delta_h=Q_h-Q_{h/r},\qquad D_h=|\delta_h|,\qquad
p_{\rm observed}=\frac{\log(D_h/D_{h/r})}{\log r}.
$$

Two levels supply one difference; at least three are needed for an order.
Four levels supply three differences and two possible order estimates.
Orders are omitted when adjacent differences change sign or reach roundoff.
Each formulation uses its own finest reference, which is a numerical approximation.
Pressure/density use conservative volume-overlap restriction; periodic Z
wrapping handles the shifted Fourier cell faces. Vorticity uses Fourier Z
and cubic X point sampling. Local errors
are volume-weighted RMS (`L2`) and maximum (`Linf`), plus relative norms with
explicit zero-reference handling. The finest-reference error is unavailable,
rather than a zero on logarithmic axes.

For the restricted reference `f_ref` and error `e=f_h-f_ref`, the local norms are

$$
L_2=\left(\frac{\sum_{\rm physical}e^2\Delta V}
                         {\sum_{\rm physical}\Delta V}\right)^{1/2},
\qquad L_\infty=\max_{\rm physical}|e|.
$$

Absolute errors retain the selected field's units. Relative errors divide by
the reference RMS or maximum norm; a zero denominator yields an unavailable
relative error. No instantaneous mean subtraction or phase alignment is
performed, and boundary-adjacent physical cells remain in the norms.

Choose common snapshot and convergence times after a pilot that identifies
growth, saturation and viscous decay, then freeze them before the final scan.
Late-time nonlinear dephasing can increase pointwise errors; every formulation
and resolution is still compared at the same declared physical time.

`--accuracy-pilot` repeats the two finest grids with tighter time/potential
solver tolerances and half the maximum timestep, then with twice the output
cadence. `accuracy-pilot.json` records observable and field changes relative to
spatial differences, with a target below 10%, and compares cumulative energy
residuals. The tightened variant divides time/recovery tolerances by ten and
halves the maximum timestep; the cadence variant doubles output intervals.
The report compares observable changes with the smallest adjacent spatial
difference and field changes with the smallest retained finest-reference L2
error. Zero spatial differences leave the corresponding ratio unavailable.
For complete energy balances, the cadence report also estimates the remaining
fine trapezoidal power-sampling error as one third of the coarse/fine residual
difference, assuming second-order sampling convergence.
Freeze the final settings after inspecting this evidence. Failure to
self-converge or a remaining energy drift is a reported numerical result.
A short smoke run does not validate temporal error for the full-duration study.

## Products

Each boundary policy produces state, energy-history and convergence figures:

- `state-*.pdf` and `state-*.png`: ion internal energy density `3*Pi/2` from
  pressure evolution (or the selected density/vorticity field) above
  perpendicular kinetic energy density from total ion energy, on each
  formulation's finest grid at four common times labelled in units of
  `1/omega_ci`. Each row has a separate colour scale shared across times,
  linear for both rows, with physical coordinates and rasterised image artists.
  Guard cells are excluded.
- `energy-*.png` (default): six panels showing thermal and perpendicular kinetic changes,
  common and native total drift, available corrected residuals and native/common discrepancy.
  Panel labels call these "shared FV energy" and "model energy"; the finest
  grid uses solid traces and coarser grids use distinct dash patterns.
- `global-convergence-*.pdf`: common `Kperp`, adjacent differences and valid
  observed orders at the specified time.
- `local-error-*.pdf`: absolute RMS and maximum error against each formulation's
  own finest reference.

`--energy-format png` saves the energy histories in one six-panel PNG per
boundary policy. `--energy-format pdf-panels` instead saves six standalone PDFs,
each with its own labels and legend:

- `energy-<boundary>-ion-thermal.pdf`
- `energy-<boundary>-electron-thermal.pdf`
- `energy-<boundary>-perpendicular-kinetic.pdf`
- `energy-<boundary>-shared-relative-change.pdf`
- `energy-<boundary>-model-relative-change.pdf`
- `energy-<boundary>-estimator-difference.pdf`

The export choice applies to the current invocation and can be changed with
`--analyse-only` without rerunning simulations. State and convergence exports
are the same in both modes; existing figure files are retained.

`energy-table` and `convergence-table` are printed and saved as CSV and JSON.
`--latex` adds tables from the same numeric data. Field units are retained in
the convergence data. Each run also retains `inventory.nc` and `analysis.json`,
including diagnostic comparisons and boundary/viscous audit information.
Inventory time is in seconds and energy-history axes use milliseconds.
CLI times use normalised `t*omega_ci`; state-image labels display the
corresponding physical time as `t=N/omega_ci`.
The energy table records initial/final native energy, signed drift, maximum
relative drift, common-estimator drift, native/common discrepancy and available
integrated outward energy and corrected-residual data. `analysis.json`
additionally reports drift relative to the peak-to-peak common perpendicular
kinetic energy when that signal is nonzero. The convergence table records physical spacing,
`Q_h`, signed adjacent and finest-reference differences, order validity and
absolute/relative local errors. Unavailable values remain explicit in the
reports and printed tables.
`figure-config.json` records plotting settings; `--style FILE` accepts a JSON
mapping of Matplotlib rcParams. Reanalysis reuses the saved rcParams; an explicit
`--style` overrides them. Fonts are embedded in the PDF and traces/labels remain
vector graphics.
