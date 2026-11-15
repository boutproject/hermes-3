# Continuum ion energy MMS

Two-ion expression sources evolve density, momentum, electron pressure,
vorticity and total ion energy through t=0.1. Default runs baseline/curvature
at 8/16 X/Z cells; long runs all four cases at 8/16/32 on two MPI ranks.
Y resolution is twice X/Z. Baseline/curvature have nonzero potential and
varying composition; skew/viscosity have zero potential and high psi harmonics.
`PROFILES` and `--variants` select coverage. The chi solve excludes the viscous
rate, and polarisation energy transport carries the viscous drift work. The
zero-potential viscous case retains vorticity diffusion and electron heating;
its analytic chi and polarisation current vanish. No separate viscous energy
flux is subtracted. Its continuum sources use this closure. Order
gates are 1.7 for 8-to-16 refinement and 1.8 for 16-to-32;
recovery/diagnostics use `1e-10`.

Every case saves 20 intervals of 0.005, retaining the end time 0.1. Errors
and orders use only the final saved output for evolved and recovered fields.
The default 300 s stall window measures lack of saved-time progress; it can therefore
detect a stopped run without waiting for its complete MMS integration. Case
provenance includes the observed output wall timestamps and their maximum
interval, sampled by the monitor at up to 1 s resolution.

`stall_seconds_by_size` selects a 3300 s stall window for 32 cells. The
8/16-cell window is 300 s; each case records its selected window and observed
output intervals. The default and long CTest time limits are 6000 and 36000 s,
with per-case limits of 3600 and 10800 s, respectively.

```sh
python tests/integrated/ion-energy-mms/runtest --exe build/hermes-3 --timeout 3600
python tests/integrated/ion-energy-mms/runtest --exe build/hermes-3 --profile long --nproc 2 --timeout 10800
python tests/integrated/ion-energy-mms/runtest --exe build/hermes-3 --variants curvature
python tests/integrated/ion-energy-mms/operator_convergence.py --output /tmp/operator-orders.json
python tests/integrated/ion-energy-mms/mms.py --all --check --check-files
```

The standalone diagnostic needs SymPy, not Hermes. `--y`/`--time` select the
analytic X-Z plane; `--run PATH` checks C++ Reynolds diagnostics. It reports
16/32/64 errors and maxima for full-vorticity MC advection and the central
Reynolds remainder. Zero-phi cases have zero analytic Reynolds terms and
undefined orders.
Runner results save errors, effective limits and per-variant times; all
failing gates stay enabled.

## Convergence checks

The production Reynolds current is `ω vE + R`: MC reconstruction for full
vorticity and a central metric remainder, retaining physical wall currents.
Polarisation work uses `φ + ψ`; wall transport uses the product of the mean
particle partition and separate means of face temperature and velocity.
The kinetic weight squares the mean velocity. `DivEExB`
advects total `Eion`. The standalone diagnostic reconstructs its perpendicular
contribution and reports the face-work diagnostic; these are separate from the
final-time MMS order gates.

MC limiting near extrema and the coupled pressure/potential recovery can lower
observed orders on these grids. Passing an integral energy budget does not imply
that the spatial or evolved-field order gates pass. Each runner invocation
records its measured errors, orders and failing checks in `result.json`.

The enabled Python `test_reynolds_wall_curvature_converges` regression checks
boundary, outer-three-cell, all-cell and interior errors independently of these
MMS runs. Its Reynolds boundary error remains non-convergent with the current
potential wall conditions and central remainder. The nonzero wall curvature
state and its failing accuracy gate remain enabled. These MMS profiles do not
establish radial wall accuracy.
