# Continuum ion energy MMS

Two-ion expression sources evolve density, momentum, electron pressure,
vorticity and total ion energy through t=0.1. Default runs baseline/curvature
at 8/16 X/Z cells; long runs all four variants at 8/16/32 on two MPI ranks.
Y resolution is twice X/Z. Baseline/curvature have nonzero potential and
varying composition; skew/viscosity have zero potential and high psi harmonics.
`PROFILES` and `--variants` select coverage. Physics, grids, end time, error
caps and order gates (1.7 / 1.8) are unchanged; recovery/diagnostics use `1e-10`.

Every case saves 20 intervals of 0.005, retaining the end time 0.1. Errors
and orders use only the final saved output for evolved and recovered fields.
The default 300 s stall window measures lack of saved-time progress; it can therefore
detect a stopped run without waiting for its complete MMS integration. Case
provenance includes the observed output wall timestamps and their maximum
interval, sampled by the monitor at up to 1 s resolution.

A standalone 32-cell baseline completed with a maximum interval of 618.45 s,
so 300 s cannot cover this grid. `stall_seconds_by_size` selects 3300 s for
32 cells only. The allowance is twice the projected curvature interval,
using the measured 16-cell curvature/baseline ratio `142.63 / 55.24`, rounded
up to 300 s. The 8/16-cell window remains 300 s; each case records its selected
window. Earlier 300 s failures are retained in the external rerun report.

At 8 cells, baseline final-time L2 errors with 20 outputs match the original
single-interval cadence for every evolved and recovered reference field.
The largest relative difference is `6.02e-13`, below the `1e-6` invariance
limit. The extra outputs provide progress without changing the solution.

```sh
python tests/integrated/ion-energy-mms/runtest --exe build/hermes-3 --timeout 3600
python tests/integrated/ion-energy-mms/runtest --exe build/hermes-3 --profile long --nproc 2 --timeout 10800
python tests/integrated/ion-energy-mms/runtest --exe build/hermes-3 --variants curvature
python tests/integrated/ion-energy-mms/operator_convergence.py --output /tmp/operator-orders.json
python tests/integrated/ion-energy-mms/mms.py --all --check --check-files
```

The standalone diagnostic needs SymPy, not Hermes. `--y`/`--time` select the
analytic X-Z plane; `--run PATH` checks C++ Reynolds diagnostics. It reports
16/32/64 term errors and maxima for both reconstructions. Zero-phi variants
have zero analytic Reynolds terms and undefined orders. Runner results save
errors, effective limits and per-variant times; all failing gates stay enabled.

## Known failures

“Before” is the split operator (`82ab0ec`); curvature has no split reference.
L2 norms exclude guards; interior norms also omit two X cells at either wall.
“At walls” means wall-band RMS is at least twice the interior RMS.

Central interior terms converge at second order; XPPM omega-vE clipping gives
orders near 1.5 and weighted currents near 1.64. Further evolved degradation
involves coupling/time evolution. Pinning adds wall-local work/weighted-curl
truncation; the annulus wall-heating integral converges faster than first order.
Baseline fine Ni crosses 1.8 from unpinned 1.80165 to pinned 1.79943, and NVj
from 1.80463 to 1.79169, through weighted wall transport/coupling.
Skew/viscosity analytic Reynolds terms vanish; coarse recovery errors include
cos(6a) and pass on refinement. Pe failures
involve unchanged shared limited transport. The new central periodic seam
correction restores roundoff work at 32/64 and preserves default XPPM behavior.
Detailed attribution accompanies each failure and the spatial table below.

### default

| Variant / variable | Before (`82ab0ec`) | After / limit | L2 all / interior | At walls? | Attribution |
| --- | --- | --- | --- | --- | --- |
| baseline / Ni (8→16) | 1.95503 | 1.66505 / 1.7 | 2.74217e-05 / 2.56975e-05 | no | Limited weighted Reynolds current; density/momentum recovery coupling |
| baseline / Pe (8→16) | 1.67305 | 1.67318 / 1.7 | 4.59725e-05 / 4.40336e-05 | no | Unchanged shared limited transport; coupled fields |
| baseline / Vort (8→16) | 1.72006 | 1.69031 / 1.7 | 8.30964e-05 / 9.52794e-05 | no | Interior omega-vE MC clipping; coupled time evolution |
| baseline / Eion (8→16) | 2.01354 | 1.62018 / 1.7 | 8.7433e-05 / 9.50571e-05 | no | Limited weighted Reynolds current; energy/recovery coupling |
| curvature / Ni (8→16) | — | 1.69440 / 1.7 | 3.00246e-05 / 2.79507e-05 | no | Limited weighted Reynolds current; density/momentum recovery coupling |
| curvature / Pe (8→16) | — | 1.61449 / 1.7 | 6.19868e-05 / 6.85143e-05 | no | Unchanged shared limited transport; coupled fields |
| curvature / Vort (8→16) | — | 1.69122 / 1.7 | 8.31519e-05 / 9.53425e-05 | no | Interior omega-vE MC clipping; coupled time evolution |
| curvature / Eion (8→16) | — | 1.66672 / 1.7 | 9.23502e-05 / 0.000100266 | no | Limited weighted Reynolds current; energy/recovery coupling |

Measured runner time: 2244.6 s; 8 checks fail.
Recovery, error caps and diagnostics on/off pass.

| Variant | Sizes | Execution time (s) |
| --- | --- | ---: |
| baseline | 8/16 | 931.5 |
| curvature | 8/16 | 1264.7 |

The diagnostics-off baseline check adds 43.3 s.

### long

| Variant / variable | Before (`82ab0ec`) | After / limit | L2 all / interior | At walls? | Attribution |
| --- | --- | --- | --- | --- | --- |
| baseline / Ni (8→16) | 1.95503 | 1.66505 / 1.7 | 2.74217e-05 / 2.56975e-05 | no | Limited weighted Reynolds current; density/momentum recovery coupling |
| baseline / Pe (8→16) | 1.67305 | 1.67318 / 1.7 | 4.59725e-05 / 4.40336e-05 | no | Unchanged shared limited transport; coupled fields |
| baseline / Vort (8→16) | 1.72006 | 1.69031 / 1.7 | 8.30964e-05 / 9.52794e-05 | no | Interior omega-vE MC clipping; coupled time evolution |
| baseline / Eion (8→16) | 2.01354 | 1.62018 / 1.7 | 8.7433e-05 / 9.50571e-05 | no | Limited weighted Reynolds current; energy/recovery coupling |
| baseline / Ni (16→32) | 1.89755 | 1.79943 / 1.8 | 7.87792e-06 / 7.64749e-06 | no | Weighted wall/coupled transport; unpinned order 1.80165 |
| baseline / NVj (16→32) | 1.80898 | 1.79169 / 1.8 | 1.82776e-06 / 1.83078e-06 | no | Weighted wall/coupled transport; unpinned order 1.80463 |
| baseline / Pe (16→32) | 1.72263 | 1.72286 / 1.8 | 1.39272e-05 / 1.36575e-05 | no | Unchanged shared limited transport; coupled fields |
| baseline / Vort (16→32) | 1.50015 | 1.03969 / 1.8 | 4.04206e-05 / 4.30931e-05 | no | Interior omega-vE MC clipping; coupled time evolution |
| baseline / Eion (16→32) | 2.29746 | 1.67095 / 1.8 | 2.74581e-05 / 2.87881e-05 | no | Limited weighted Reynolds current; energy/recovery coupling |
| skew / Ni (8→16) | 1.76549 | 1.68311 / 1.7 | 2.5033e-05 / 2.22789e-05 | no | Recovery/limited coupled transport; analytic Reynolds terms vanish |
| skew / Pe (8→16) | 1.68983 | 1.68984 / 1.7 | 3.99162e-05 / 3.56076e-05 | no | Unchanged shared limited transport; coupled fields |
| skew / psi (8→16) | 1.38208 | 1.36537 / 1.7 | 0.000160965 / 0.000184545 | no | High psi harmonics; pre-asymptotic recovery, passes 16→32 |
| skew / Tion (8→16) | 1.37634 | 1.37209 / 1.7 | 0.00013687 / 0.000157962 | no | High psi harmonics; pre-asymptotic recovery, passes 16→32 |
| skew / Pi (8→16) | 1.38164 | 1.36508 / 1.7 | 0.000169403 / 0.000194325 | no | High psi harmonics; pre-asymptotic recovery, passes 16→32 |
| skew / pi (8→16) | 1.38442 | 1.36715 / 1.7 | 0.000135758 / 0.000155278 | no | High psi harmonics; pre-asymptotic recovery, passes 16→32 |
| skew / pj (8→16) | 1.38810 | 1.37441 / 1.7 | 3.43345e-05 / 3.95173e-05 | no | High psi harmonics; pre-asymptotic recovery, passes 16→32 |
| viscosity / Pe (8→16) | 1.68991 | 1.68991 / 1.7 | 3.98966e-05 / 3.55795e-05 | no | Unchanged shared limited transport; coupled fields |
| viscosity / Vort (8→16) | 1.65083 | 1.64784 / 1.7 | 1.81269e-05 / 2.04523e-05 | no | High psi harmonics; pre-asymptotic recovery, passes 16→32 |
| viscosity / psi (8→16) | 1.33499 | 1.31792 / 1.7 | 0.00015941 / 0.000182715 | no | High psi harmonics; pre-asymptotic recovery, passes 16→32 |
| viscosity / Tion (8→16) | 1.33086 | 1.31966 / 1.7 | 0.0001343 / 0.000154992 | no | High psi harmonics; pre-asymptotic recovery, passes 16→32 |
| viscosity / Pi (8→16) | 1.33452 | 1.31743 / 1.7 | 0.000167734 / 0.000192364 | no | High psi harmonics; pre-asymptotic recovery, passes 16→32 |
| viscosity / pi (8→16) | 1.33759 | 1.32062 / 1.7 | 0.000134546 / 0.000153845 | no | High psi harmonics; pre-asymptotic recovery, passes 16→32 |
| viscosity / pj (8→16) | 1.34338 | 1.32631 / 1.7 | 3.38903e-05 / 3.90005e-05 | no | High psi harmonics; pre-asymptotic recovery, passes 16→32 |
| curvature / Ni (8→16) | — | 1.69440 / 1.7 | 3.00246e-05 / 2.79507e-05 | no | Limited weighted Reynolds current; density/momentum recovery coupling |
| curvature / Pe (8→16) | — | 1.61449 / 1.7 | 6.19868e-05 / 6.85143e-05 | no | Unchanged shared limited transport; coupled fields |
| curvature / Vort (8→16) | — | 1.69122 / 1.7 | 8.31519e-05 / 9.53425e-05 | no | Interior omega-vE MC clipping; coupled time evolution |
| curvature / Eion (8→16) | — | 1.66672 / 1.7 | 9.23502e-05 / 0.000100266 | no | Limited weighted Reynolds current; energy/recovery coupling |
| curvature / Pe (16→32) | — | 1.59572 / 1.8 | 2.05088e-05 / 2.15733e-05 | no | Unchanged shared limited transport; coupled fields |
| curvature / Vort (16→32) | — | 1.04073 / 1.8 | 4.04187e-05 / 4.30908e-05 | no | Interior omega-vE MC clipping; coupled time evolution |
| curvature / Eion (16→32) | — | 1.71262 / 1.8 | 2.81766e-05 / 2.95534e-05 | no | Limited weighted Reynolds current; energy/recovery coupling |

Measured runner time: 8405.8 s; 30 checks fail.
Recovery, error caps and diagnostics on/off pass.

| Variant | Sizes | Execution time (s) |
| --- | --- | ---: |
| baseline | 8/16/32 | 2719.1 |
| skew | 8/16/32 | 240.8 |
| viscosity | 8/16/32 | 761.9 |
| curvature | 8/16/32 | 4647.5 |

The diagnostics-off baseline check adds 23.8 s.

### Spatial operator attribution

Orders below are 32→64 on the analytic Y=2.91, t=0.1 plane. Baseline and curvature agree. Full term errors and maxima for 16/32/64 are saved by the diagnostic.

| Term | Central all / interior | XPPM all / interior | Finest XPPM maximum |
| --- | --- | --- | --- |
| Lambda | 2.00 / 2.05 | 2.00 / 2.05 | (4.92, 0.08), interior |
| omega_vE | 1.99 / 2.04 | 1.46 / 1.51 | (4.92, 2.42), interior |
| Q | 1.97 / 2.02 | 1.97 / 2.02 | (3.83, 7.58), interior |
| Lambda_U | 1.98 / 2.03 | 1.80 / 1.85 | (4.92, 4.30), interior |
| DivJReynolds | 1.98 / 2.02 | 1.59 / 1.64 | (4.92, 2.42), interior |
| DivEReynolds | 1.96 / 2.01 | 1.96 / 2.01 | (5.70, 2.11), interior |
| work | 0.98 / 2.02 | 0.93 / 2.02 | (0.08, 9.14), wall band |
| transport_Ni | 1.97 / 2.02 | 1.59 / 1.64 | (5.08, 7.42), interior |
| S_grad_Ni | -0.53 / 2.03 | -0.53 / 2.03 | (0.08, 7.58), wall band |
| transport_Nj | 1.90 / 2.02 | 1.56 / 1.64 | (4.92, 2.42), interior |
| S_grad_Nj | -0.53 / 2.03 | -0.53 / 2.03 | (0.08, 7.58), wall band |
| transport_NVi | 0.56 / 2.02 | 0.71 / 1.64 | (9.92, 7.58), wall band |
| S_grad_NVi | -0.54 / 2.03 | -0.54 / 2.03 | (9.92, 2.58), wall band |
| transport_NVj | 1.49 / 2.02 | 1.38 / 1.64 | (4.92, 2.42), interior |
| S_grad_NVj | -0.54 / 2.03 | -0.54 / 2.03 | (9.92, 4.92), wall band |
| transport_E_adv | 1.93 / 2.02 | 1.57 / 1.64 | (4.92, 2.27), interior |
| S_grad_E_adv | -0.43 / 2.03 | -0.43 / 2.03 | (9.92, 2.42), wall band |

With the finer output cadence, the default validation took 2544.14 s: all
five runs completed, recovery and error caps passed, and the same eight order
gates failed. The longest output interval in that validation was 162.1 s.
Default CTest timeout is 6000 s, exceeding twice that profile time; long is
36000 s. Per-case limits remain 3600 / 10800 s.

The four C++ Reynolds diagnostics agree with the mirror within `4.2e-17` on the baseline-8 output. Timing includes concurrent validation load.
