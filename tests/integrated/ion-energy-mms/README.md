# Continuum ion energy MMS

All ion-energy MMS variants use this directory and the same runner. Select one
with `--variant`. The default is `parallel_current`. The runner uses committed expressions and
does not require SymPy; generating and symbolically checking inputs does.

| Variant | Purpose | Acceptance region |
| --- | --- | --- |
| `parallel_current` | 3D independent electron momentum and parallel current; φ and ψ have zero radial second derivatives at both walls | Full physical domain |
| `parallel_current_wall` | Same prescribed current, densities, velocities and electron temperature; ψ has nonzero radial curvature at both walls | Fixed interior: normalized X ∈ [0.25, 0.75] |
| `baseline` | Zero-current two-ion slab with potential and Reynolds stress | Full physical domain |
| `skew` | Zero-current skew-metric slab with zero potential | Full physical domain |
| `viscosity` | Zero-current viscous slab with zero potential | Full physical domain |
| `curvature` | Zero-current slab with magnetic curvature | Full physical domain |

The first two cases isolate the radial-wall Reynolds error. The pressure
potential uses P(a)=(1−cos(2a))² in the smooth case and P(a)=1−cos(2a) in the
curved case, where a=πX/L. Analytic chi inverses are supplied for both. Both envelope
Z-dependent ion density, velocity and temperature harmonics with sin²(πX/L),
so their wall values are independent of Z. This makes the continuum solution
compatible with the separate Z averaging of wall transport weights. Their L2 and L∞
order limits are **1.4 and 0.9**, respectively, consistent with MC limiting near
extrema. The curved-wall variant checks these orders in the fixed central half
of the radial domain, and reports global errors and orders without gating them.
A fixed physical interval avoids changing the sampled region with refinement.
The smooth-wall variant checks global orders. The condition concerns the
Reynolds potentials, not the curvature of every primitive field.

L2 is a volume-weighted RMS; L∞ is the maximum absolute error. Both exclude
mesh guard cells. Errors are measured at the final saved time for every evolved
and recovered field, including electron momentum and chi. Large chi errors are
diagnostic-only; finite errors and convergence remain required. Chi reference
RMS and relative RMS error are saved. Other L2 error caps, source ownership,
energy recovery, current continuity and nonzero-current checks remain enabled.
The four zero-current variants retain L2 order limits of 1.7 on the first
refinement and 1.8 thereafter. These historical gates remain active, including
known convergence failures; consolidation does not establish that those tests
pass.

## Running and publication output

```sh
# Two levels (8 and 16 X/Z cells), one variant per invocation
./tests/integrated/ion-energy-mms/runtest --variant parallel_current
./tests/integrated/ion-energy-mms/runtest --variant parallel_current_wall

# Publication figures: four levels by default (8, 16, 32, 64)
./tests/integrated/ion-energy-mms/runtest --variant parallel_current --plot --nproc 2 --timeout 10800
./tests/integrated/ion-energy-mms/runtest --variant parallel_current_wall --plot --nproc 2 --timeout 10800

# An explicit number of levels overrides the plotting default
./tests/integrated/ion-energy-mms/runtest --variant parallel_current_wall --plot --levels 3
./tests/integrated/ion-energy-mms/runtest --variant baseline --levels 3 --verify

# Symbolic identities, input freshness, and independent operator diagnostics
python tests/integrated/ion-energy-mms/mms.py --all --check --check-files
python tests/integrated/ion-energy-mms/operator_convergence.py --output /tmp/operator-orders.json
```

Use `--exe` to select a build and `--output-root` to keep results in a chosen
parent directory. Each level doubles X, Y and Z resolution, starting at 8 × 16
× 8 physical cells. At least two levels are required. Parallel-current cases
end at t=0.02 with four output intervals; legacy cases end at t=0.1 with twenty.
Both use output intervals of 0.005. Independent electron momentum introduces a
fast timescale, so fine publication runs can take hours. `--timeout` is per
solver launch (default 3600 s). Lack of saved-time progress is checked over 600 s
at 8 cells, 1800 s at 16 cells and 3300 s at 32 cells and above.

CTest registers both parallel-current variants with `--verify` on two ranks.
Their three-level registrations carry the `long` label, alongside separate
three-level registrations for each legacy variant. Plotting remains opt-in.

`--plot` needs matplotlib and writes a CSV table plus per-field PDF and SVG
figures with domain/interior L2 and L∞ curves and order-one/order-two reference
lines. SVG text remains editable and PDF fonts are embedded. The loader,
reconstructed quantities and volume integration use the shared xhermes-backed
`tools/energy_analysis` adapters; plotting consumes the recorded norms.

Without plotting, a failed gate stops before the next launch. Plot runs retain
spatial accuracy failures and collect the remaining requested resolutions;
execution or physics failures still stop immediately. Figures from completed
levels and `result.json` are saved on either path, and any failed check gives a
nonzero exit status. This allows publication of degraded global convergence
without suppressing failures of the selected acceptance region.

`--verify` additionally checks tighter temporal integration at the second
resolution and serial/MPI agreement at the coarsest resolution for the parallel
cases. For a zero-current variant it checks diagnostics-on/off agreement.
Dependent checks stop after a failure; verification is skipped if a plot run
has failed spatial acceptance. Continuum forcing is undamped and does not
cancel numerical `phi_dissipation`. Damping is enabled in both parallel inputs.

`--reuse-root PATH` accepts completed variant-resolution runs only when input,
executable and override hashes match. `--internal-progress` enables optional
CVODE timestep hooks in a diagnostic build. Neither is needed for normal runs.

The Reynolds current is `ω vE + R`, with MC advection and a central metric
remainder. The independent operator regressions demonstrate bounded first-row
error for curved walls and the absence of that leading error for zero-curvature
walls. Coupled MMS curves must still be measured; conservation and symbolic
identities alone do not establish final-time convergence.
