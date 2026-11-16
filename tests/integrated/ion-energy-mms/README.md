# Continuum ion energy MMS

The ion-energy component pairs explicitly with conservative `polarisation_drift`.
Auxiliary AC boundary settings are owned by the drift section. These existing
variants retain their zero-current electron closure.

The additional generated `parallel_current` variant evolves independent
electron momentum. Its three-level, temporal-refinement and MPI accuracy gates
are registered separately as `ion-energy-parallel-mms`; see the
[parallel MMS README](../ion-energy-parallel-mms/README.md). Existing default
and long profile selections and their thresholds remain unchanged.

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
recovery/diagnostics use `1e-10`. Field-error limits are separate entries in
`criteria.json["errors"]`, e.g. Vort `8.165673925347495e-4` and Eion
`1.1908840233405635e-3`.

Every case saves 20 intervals of 0.005, retaining the end time 0.1. Errors
and orders use only the final saved output for evolved and recovered fields.
The current pre-extraction baseline fails the default order gates for Ni, Pe,
Vort and Eion in both baseline and curvature variants, with measured orders
between 1.375 and 1.694 against 1.7. Recovery and absolute-error gates pass.
Comparison with `criteria.json["observed"]` shows all eight historical orders
were already below 1.7. Six fresh values reproduce them within 0.001; baseline
and curvature vorticity decline from about 1.69 to 1.38 before extraction and
need a separate investigation. The older passing-suite narrative conflicts
with these recorded criteria. The task's named deferral covers only these
eight default 8-to-16 order gates, retaining raw failures and all thresholds.
Long-profile and new parallel-current MMS
accuracy requirements remain in force.
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

The Python `test_reynolds_curved_wall_limits_global_convergence` regression
retains the nonzero wall-curvature state and current, checking first-order
domain mean absolute error, near-second-order interior RMS and a bounded
nonvanishing wall error. `test_reynolds_zero_wall_curvature_converges` checks
second-order wall maximum and domain mean absolute accuracy with production
reflection and MC donors. Its centred-donor control also checks second-order
domain RMS/maximum accuracy. These operator checks are independent of the
coupled MMS gates; passing an integral budget alone does not establish wall
accuracy.
