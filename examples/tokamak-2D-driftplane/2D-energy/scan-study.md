# Interchange resolution study

`run_scan.py` runs matched pressure and total-ion-energy formulations
and retains the inputs, time-zero restart files, raw dumps, logs and scalar
analysis. Energy drift and spatial self-convergence are separate measurements:
conserving total energy does not establish spatial accuracy. Routine conservation
checks use `tests/integrated/ion-energy-closed`; the optional 3D convergence
study uses `tests/integrated/ion-energy-mms`.

Install the packages in `tools/energy_analysis/requirements.txt` plus
`matplotlib` into your Python environment:

```sh
python -m pip install -r tools/energy_analysis/requirements.txt matplotlib
```

Run the study from the repository root:

```sh
python examples/tokamak-2D-driftplane/2D-energy/run_scan.py \
  --exe build/hermes-3 --phi-boundary both --accuracy-pilot \
  --output-root /path/to/results

# Bounded validation: 8² and 16², fixed and relaxing boundaries, viscosity=0.
python examples/tokamak-2D-driftplane/2D-energy/run_scan.py --smoke

# Recreate plots and tables without executing Hermes.
python examples/tokamak-2D-driftplane/2D-energy/run_scan.py \
  --analyse-only /path/to/results/interchange-TIMESTAMP --latex

# Export each energy-history panel as a standalone PDF.
python examples/tokamak-2D-driftplane/2D-energy/run_scan.py \
  --analyse-only /path/to/results/interchange-TIMESTAMP --energy-format pdf-panels
```

## Run configuration

The default matrix has eight runs: both formulations on four physical grids
(32², 64², 128², 256²), with fixed phi boundaries. `--phi-boundary both`
produces sixteen runs. Cases run sequentially and each owns its inputs and
outputs. The full matrix is an explicit longer study; `--smoke` is
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
| `--custom-D` | `1.19` | Positive normalised constant particle and thermal-energy diffusivity, common to both formulations |
| `--end-time` | `20000` | Duration in normalised time, `t*omega_ci` |
| `--outputs` | `200` | Regular output intervals, plus the initial record; at least four |
| `--snapshot-times` | `0, T/4, T/2, T` | Four normalised snapshot times, where `T` is the duration |
| `--convergence-time` | `T/2` | Normalised time for global and local convergence |
| `--field` | `ion_pressure` | Local convergence field; also `density` or `vorticity` |
| `--nproc` | `1` | MPI ranks, decomposing X only |
| `--timeout` | `7200` | Wall-clock limit [s] for each startup or simulation launch |
| `--exe` | `hermes-3` or `build/hermes-3` | Repository executable, preferring the former when present |
| `--output-root` | `results` beside the script | Parent of a new `interchange-TIMESTAMP` study |

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
Laplacian boundary flags are supplied only to the pressure formulation.
Run the maintained script from the repository root.

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

Every run gets its own generated `BOUT.inp` from `scan-pressure.inp.template`
or `scan-total-energy.inp.template`. Component lists, evolving-field sections
and operator settings are explicit in each template; the runner substitutes
scan parameters and the selected wall condition. The manifest hashes both
templates when checking resume compatibility.
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
missing records and irregular output axes in completed runs. It skips cases
without completed outputs and analyses the available subset, including a
single formulation or grid. Every retained input must still match its manifest
hash. No completed runs is an error. Reanalysing an older study preserves its
original physics; adding diffusion requires a new study because configuration
and template hashes prevent resuming old inputs with changed physics.

A nonzero process exit can occur after the final output, for example during
MPI finalisation. The runner checks every rank's complete expected time axis
and finite required physical fields in the final record. If these checks pass,
it accepts the outputs and records `outputs_complete=true` in `status.json`,
preserving the nonzero exit code and printing a warning. Startup recovery also
requires readable restart files for every rank. Resume and reanalysis perform
the same check for previously failed runs. Incomplete or unreadable outputs
remain failures; a keyboard interrupt stops the current invocation even when
its outputs are complete.

## Extending and selecting a study

To add resolution levels without repeating all the original simulation flags,
use `--extend`. The number of levels is the **new total**, not the number to
add. All saved simulation settings and the original executable path are
inherited; only an increased level count is accepted.

```sh
python examples/tokamak-2D-driftplane/2D-energy/run_scan.py \
  --extend /path/to/results/interchange-TIMESTAMP --levels 5 \
  --allow-runner-change
```

The manifest retains existing cases and records the appended grids. Completed
cases are reused; remaining cases run normally. Once the whole matrix is
complete, figures and tables are regenerated, using each formulation's new
finest grid as its convergence reference. Decreasing the level count or
changing physical settings is rejected. The runner-version override is needed
only when the script has changed since the study was created or last resumed.

`--cases CASE [CASE ...]` selects particular directory names for local
execution, such as `fixed-pressure-128x128` and `fixed-total_energy-128x128`.
Completed cases are skipped. When cases remain pending, use `--resume`
to complete the matrix or `--analyse-only` to plot the completed subset.

## Physical setup and startup

Both models evolve hydrogen density, impose quasineutral electrons, and evolve
electron pressure. The pressure formulation also evolves ion pressure and uses
full Reynolds stress and legacy polarisation advection. The total-energy model
evolves `Eion` and recovers the common ion temperature. Both include
`classical_diffusion` in its positive `custom_D` mode, diffusing particle density
and thermal energy with the same constant coefficient. For each species, the
thermal-energy source is `1.5*D*laplacian(P)`. Physical diffusion faces are closed;
Z is periodic. The coefficient stays fixed across grids. The normalisation is
`D_SI = custom_D * rho_s0^2 * Omega_ci`; here `rho_s0^2 * Omega_ci = 14.2857 m²/s`,
so the default `1.19` gives 17 m²/s. A positive override avoids requiring a
collision-frequency component in this setup.
There is no parallel flow, additional collisional conduction, prescribed source
or hyperdiffusion. Both use reference atomic mass 1.
The box is 5 cm × 5 cm, with R=1.5 m, B=0.35 T, Nnorm=2e18 m⁻³ and Tnorm=5 eV.
Both templates supply consistent SI metric pairs explicitly and use
`normalise_metric=true`; the legacy `recalculate_metric` path mixes B units
when constructing its Z metric for this setup. Guard widths stay fixed when
physical grid counts are refined: `MXG=2`, `MYG=0`, `mesh:ny=1` and
`mesh:nx=nx_physical+4`. Thus the default X mesh sizes are 36, 68, 132 and 260.
X has physical walls, Z is periodic, and the constant curvature is
`bxcvz=1/Rxy^2`. Cell widths are recomputed at each resolution so the box and
perturbation wavelengths remain fixed.

Initial density is uniform and both species have
`P=1+0.1*S(x)*sin(2*pi*x)*sin(2*pi*x-z)`, with `x` between 0 and 1. The
envelope `S` is zero within 0.1 of either radial wall, one at distances of
0.2 or greater, and uses the quintic `s^3*(10-15*s+6*s^2)` between them,
where `s=(distance-0.1)/0.1`. Its first and second derivatives vanish at
both taper endpoints. The buffer and taper widths are fractions of the box,
fixed across resolutions; the central perturbation is unchanged. The
templates and restart builder use the same profile.
The shared physical potential is sampled from the analytic continuum ion-only
zero-vorticity solution: the harmonic extension is `Phi=1`, so `phi=1-P`.
Both models have zero phi midpoints. The pressure formulation uses homogeneous Neumann
ion-pressure guards; total ion energy version 1 uses Neumann guards on ion
pressure and its pressure potential, local Φ wall work and local transport
weights. Electron pressure retains
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
cell-centred heating defect is reconstructed and reported. The constant
wall buffers satisfy both initial pressure boundary conditions and give
zero initial wall gradients of phi and the pressure potential. The legacy
electron mass contribution leaves a small, smooth bulk vorticity; total ion
energy has zero initial vorticity and perpendicular kinetic energy to roundoff.
The previous profile vanished at the walls but had a nonzero radial pressure
derivative. Combining Neumann pressure guards with zero-midpoint phi guards
gave the legacy generalised potential a grid-scale kink, producing a
first-cell vorticity proportional to inverse radial spacing. Increasing
diffusion suppressed this sheet but also damped the bulk dynamics.
Boundary-adjacent physical cells are retained in all primary error norms.
These differences are recorded in every manifest.

Use `--custom-D 0.001` (0.0143 m²/s) to study weakly diffusive dynamics;
the historical default `1.19` is strongly diffusive. A fresh study is required
for the buffered initial profile because the templates and restart fields have
changed. Existing data can still be inspected with `--analyse-only`.

The local runner limits each launch to `--timeout` seconds and records
`termination_reason=timeout` in `status.json` when this limit is reached.
A timeout does not by itself indicate a numerical exception; inspect
`run.log` and the last saved time before choosing a longer limit.

Weak-diffusion validation covered both fixed-boundary formulations to
`t=10000` on 32²/64² and to `t=100` on 128²/256², with viscosity 1 m²/s.
Full-duration 256² convergence remains unverified. Later wall vorticity
can develop in both models; removing the initial sheet does not ensure
that evolving wall layers converge.

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
Classical diffusion closes its physical faces and contributes no outward power.
For diffusion runs, the full local RHS diagnostic reconstruction is unavailable
because it does not yet include interior diffusion sources; `analysis.json`
records this limitation when `--diagnose` is enabled. The global wall-power
audit remains available. Older studies without diffusion retain the full
operator cross-checks.
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
Extending a study recomputes these differences over all completed grids and
uses the new finest grid for local errors. Five levels provide three possible
order estimates, regardless of the original level count.
Orders use the absolute adjacent differences even when their signs change.
`order_status` flags sign changes as `oscillatory differences`; a positive order
means their magnitudes shrink under refinement. This measures the decay of
adjacent differences without assuming a single signed leading error term.
Orders are omitted when either adjacent difference reaches roundoff.
Completed subsets with unequal refinement ratios also omit adjacent orders;
equally spaced refinement triples use their actual ratio.
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

The thermal-gain plot uses the full duration, independently of
`--convergence-time`. Its signed observable is

$$
\Delta U_i(h)=U_i(T;h)-U_i(0;h),\qquad
U_i=\sum_{\rm physical}\frac32 P_i\,\Delta V.
$$

Pressure and total-ion-energy gains are plotted together when available, in
joules against `h=Delta X=Lrad/nx_physical` on a logarithmic axis; Z refines in
the same proportion. The ordinate remains linear to retain negative gains.
For at least three distinct ion-energy grids, the script fits

$$
\Delta U_i(h)=\Delta U_{i,0}+A(h/h_{\rm ref})^p,\qquad p>0,
$$

where `h_ref` is the coarsest available spacing. The fitted offset estimates the
zero-spacing thermal gain; it is not assumed to vanish. The dashed curve extends
below the smallest measured spacing and the dotted asymptote labels the limit,
since zero itself cannot appear on a logarithmic axis. The legend gives the
convergence order. The JSON records the offset, amplitude, order, reference
spacing, RMS fit residual and R². Orders are searched over `[1e-4,16]`; boundary
optima, oscillatory gains and roundoff-dominated sequences produce an explicit
unavailable fit. Three points determine three parameters without testing the
power-law model; more grids and the accuracy pilot are needed to assess the
observed convergence regime. The two-grid smoke profile produces the PDF with
an unavailable fit annotation.

## Products

Each boundary policy produces state, energy-history and convergence figures:

- `state-*.pdf` and `state-*.png`: ion internal energy density `3*Pi/2` from
  pressure evolution (or the selected density/vorticity field) above
  perpendicular kinetic energy density from total ion energy, on each
  formulation's finest grid at four common times labelled in units of
  `1/omega_ci`. Each row has a separate colour scale shared across times,
  linear for both rows, with physical coordinates and rasterised image artists.
  Guard cells are excluded. Only available formulation rows are drawn during
  partial reanalysis.
- `energy-*.png` (default): six panels showing thermal and perpendicular kinetic changes,
  common and native total drift, available corrected residuals and native/common discrepancy.
  Panel labels call these "shared FV energy" and "model energy"; the finest
  grid uses solid traces and coarser grids use distinct dash patterns.
- `global-convergence-*.pdf`: common `Kperp`, adjacent differences and valid
  observed orders at the specified time.
- `ion-thermal-gain-<boundary>.pdf` and matching `.json`: signed final-minus-initial
  ion thermal energy for available formulations, plus the ion-energy continuum
  fit and its status. Always exported, independently of `--energy-format`.
- `local-error-*.pdf`: absolute RMS and maximum error against each formulation's
  own finest reference. Each panel includes a power-law fit through the two
  highest-resolution total-ion-energy error points, with its order shown to
  two significant figures in the legend. The finest reference has no error
  point and is excluded from the fit.

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
The printed convergence table includes `order_status`, flagging oscillatory
differences and explaining unavailable orders such as roundoff or an
insufficient refinement triple.
`--latex` adds tables from the same numeric data. Field units are retained in
the convergence data. Each run also retains `inventory.nc` and `analysis.json`,
including diagnostic comparisons and boundary/viscous audit information.
Inventory time is in seconds and energy-history axes use milliseconds.
CLI times use normalised `t*omega_ci`; state-image labels display the
corresponding physical time as `t=N/omega_ci`.
The energy table records initial/final ion thermal energy and its signed gain,
initial/final native energy, signed drift, maximum
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
