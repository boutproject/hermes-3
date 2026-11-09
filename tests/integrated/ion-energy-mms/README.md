# Continuum ion energy MMS

Expression-only sources evolve two ions, zero-current electrons, vorticity
and total ion energy through time 0.1. The restored baseline has nonzero phi,
E×B transport, Reynolds stress, varying composition and finite parallel flows.
The skew metric (`g13=1/4`, `g33=17/16`) and viscosity (`nu=0.02`) variants
retain the original wall-compatible zero-phi solutions and varying composition.
All variants run at 8/16 X/Z cells with twice as many Y cells; the long profile
adds 32 and uses two MPI ranks. The new curvature variant retains nonzero phi
and composition but chooses a finite cosine mode for Pe + Pi, giving an
analytic inverse for its diamagnetic forcing. L2 orders remain 1.7 (8→16) and 1.8 (16→32).
Finest-grid limits are ten times the measured errors. Recovery and diagnostics
on/off comparisons use a fixed `1e-10` floor.

The baseline, skew and viscosity analytic fields are retained. The skew-current
identity is checked symbolically by `mms.py --check`; its vorticity divergence
agrees with the old split. Baseline sources were regenerated because the kinetic
energy flux and the weighted polarisation current change. The added curvature
variant uses `Pe = 4 + 0.02*cos(2a)*sin(theta) - Pi` and adds the corresponding
cosine eigenmode to chi. Preflight checks the full inversion and wall conditions.

```sh
python tests/integrated/ion-energy-mms/runtest --exe build/hermes-3 --timeout 3600
python tests/integrated/ion-energy-mms/runtest --exe build/hermes-3 --profile long --nproc 2 --timeout 10800
python tests/integrated/ion-energy-mms/mms.py --all --check --check-files
```

The runner needs no SymPy or plotting packages; SymPy is used for generation
and preflight. Results save all-cell L2 errors and interior L2 errors excluding
two physical cells next to each X face, plus effective limits and run times.

## Known failures

Every failing convergence check is listed below; all checks remain enabled.
Orders use all physical cells, excluding guards. Interior errors also omit two
physical cells next to each X face. “At walls” means the wall-band RMS error is
at least twice the interior RMS error. “Before” uses the isolated original
operator on the same grids; curvature is new and has no before comparison.
The PR validation description lists every field's before/after orders.

Baseline Vort convergence **worsened** from 1.72006 to 1.69031 at 8→16 and
from 1.50015 to 1.03970 at 16→32 (two-rank long runs). Its former 1.50 order
was traced to MC limiting at an interior extremum, with boundary error about
ten times below the domain maximum (`3831b4d`). That attribution does not fully
explain the new degradation. The new 32-cell error maximum remains in the
interior; the two-cell wall-band maximum is only 6.4% of the domain maximum.
Pe orders near 1.67 come from limited advection
at extrema. The coarse skew/viscosity psi, Tion and pressure failures are
pre-asymptotic: their psi includes harmonics through cos(6a), and these fields
pass at 16→32. Coarse viscosity Vort likewise requires refinement. The new
baseline Ni/Eion and skew Ni failures involve coupled transport truncation;
their individual-term causes have not been isolated. Curvature shares the
baseline potential and pressure-potential harmonics. Its fine Vort order shows
the same degradation as baseline; the individual-term causes of its Pe and
Eion failures remain unisolated. Its constructed physics and thresholds are
unchanged.

### default

| Variant / variable | Order / limit | Before | L2 all / interior | At walls? |
| --- | --- | --- | --- | --- |
| baseline / Ni (8→16) | 1.66512 / 1.7 | 1.95503 | 2.74206e-05 / 2.56974e-05 | no |
| baseline / Pe (8→16) | 1.67318 / 1.7 | 1.67305 | 4.59727e-05 / 4.40336e-05 | no |
| baseline / Vort (8→16) | 1.69031 / 1.7 | 1.72006 | 8.30964e-05 / 9.52794e-05 | no |
| baseline / Eion (8→16) | 1.62040 / 1.7 | 2.01354 | 8.74205e-05 / 9.50571e-05 | no |
| skew / Ni (8→16) | 1.68311 / 1.7 | 1.76549 | 2.5033e-05 / 2.22789e-05 | no |
| skew / Pe (8→16) | 1.68984 / 1.7 | 1.68983 | 3.99162e-05 / 3.56076e-05 | no |
| skew / psi (8→16) | 1.36537 / 1.7 | 1.38208 | 0.000160965 / 0.000184545 | no |
| skew / Tion (8→16) | 1.37209 / 1.7 | 1.37634 | 0.00013687 / 0.000157962 | no |
| skew / Pi (8→16) | 1.36508 / 1.7 | 1.38164 | 0.000169403 / 0.000194325 | no |
| skew / pi (8→16) | 1.36715 / 1.7 | 1.38442 | 0.000135758 / 0.000155278 | no |
| skew / pj (8→16) | 1.37441 / 1.7 | 1.38810 | 3.43345e-05 / 3.95173e-05 | no |
| viscosity / Pe (8→16) | 1.68991 / 1.7 | 1.68991 | 3.98966e-05 / 3.55795e-05 | no |
| viscosity / Vort (8→16) | 1.64784 / 1.7 | 1.65083 | 1.81269e-05 / 2.04523e-05 | no |
| viscosity / psi (8→16) | 1.31792 / 1.7 | 1.33499 | 0.00015941 / 0.000182715 | no |
| viscosity / Tion (8→16) | 1.31966 / 1.7 | 1.33086 | 0.0001343 / 0.000154992 | no |
| viscosity / Pi (8→16) | 1.31743 / 1.7 | 1.33452 | 0.000167734 / 0.000192364 | no |
| viscosity / pi (8→16) | 1.32062 / 1.7 | 1.33759 | 0.000134546 / 0.000153845 | no |
| viscosity / pj (8→16) | 1.32631 / 1.7 | 1.34338 | 3.38903e-05 / 3.90005e-05 | no |
| curvature / Ni (8→16) | 1.69445 / 1.7 | — | 3.00236e-05 / 2.79507e-05 | no |
| curvature / Pe (8→16) | 1.61449 / 1.7 | — | 6.19868e-05 / 6.85143e-05 | no |
| curvature / Vort (8→16) | 1.69122 / 1.7 | — | 8.31519e-05 / 9.53425e-05 | no |
| curvature / Eion (8→16) | 1.66692 / 1.7 | — | 9.23384e-05 / 0.000100266 | no |

Measured runner time: 2827.2 s; 22 order checks fail. Recovery, error caps and diagnostics on/off checks pass.

### long

| Variant / variable | Order / limit | Before | L2 all / interior | At walls? |
| --- | --- | --- | --- | --- |
| baseline / Ni (8→16) | 1.66512 / 1.7 | 1.95503 | 2.74206e-05 / 2.56974e-05 | no |
| baseline / Pe (8→16) | 1.67318 / 1.7 | 1.67305 | 4.59727e-05 / 4.40336e-05 | no |
| baseline / Pe (16→32) | 1.72286 / 1.8 | 1.72263 | 1.39274e-05 / 1.36575e-05 | no |
| baseline / Vort (8→16) | 1.69031 / 1.7 | 1.72006 | 8.30964e-05 / 9.52794e-05 | no |
| baseline / Vort (16→32) | 1.03970 / 1.8 | 1.50015 | 4.04205e-05 / 4.3093e-05 | no |
| baseline / Eion (8→16) | 1.62040 / 1.7 | 2.01354 | 8.74205e-05 / 9.50571e-05 | no |
| baseline / Eion (16→32) | 1.67564 / 1.8 | 2.29746 | 2.73651e-05 / 2.87881e-05 | no |
| skew / Ni (8→16) | 1.68311 / 1.7 | 1.76549 | 2.5033e-05 / 2.22789e-05 | no |
| skew / Pe (8→16) | 1.68984 / 1.7 | 1.68983 | 3.99162e-05 / 3.56076e-05 | no |
| skew / psi (8→16) | 1.36537 / 1.7 | 1.38208 | 0.000160965 / 0.000184545 | no |
| skew / Tion (8→16) | 1.37209 / 1.7 | 1.37634 | 0.00013687 / 0.000157962 | no |
| skew / Pi (8→16) | 1.36508 / 1.7 | 1.38164 | 0.000169403 / 0.000194325 | no |
| skew / pi (8→16) | 1.36715 / 1.7 | 1.38442 | 0.000135758 / 0.000155278 | no |
| skew / pj (8→16) | 1.37441 / 1.7 | 1.38810 | 3.43345e-05 / 3.95173e-05 | no |
| viscosity / Pe (8→16) | 1.68991 / 1.7 | 1.68991 | 3.98966e-05 / 3.55795e-05 | no |
| viscosity / Vort (8→16) | 1.64784 / 1.7 | 1.65083 | 1.81269e-05 / 2.04523e-05 | no |
| viscosity / psi (8→16) | 1.31792 / 1.7 | 1.33499 | 0.00015941 / 0.000182715 | no |
| viscosity / Tion (8→16) | 1.31966 / 1.7 | 1.33086 | 0.0001343 / 0.000154992 | no |
| viscosity / Pi (8→16) | 1.31743 / 1.7 | 1.33452 | 0.000167734 / 0.000192364 | no |
| viscosity / pi (8→16) | 1.32062 / 1.7 | 1.33759 | 0.000134546 / 0.000153845 | no |
| viscosity / pj (8→16) | 1.32631 / 1.7 | 1.34338 | 3.38903e-05 / 3.90005e-05 | no |
| curvature / Ni (8→16) | 1.69445 / 1.7 | — | 3.00236e-05 / 2.79507e-05 | no |
| curvature / Pe (8→16) | 1.61449 / 1.7 | — | 6.19868e-05 / 6.85143e-05 | no |
| curvature / Pe (16→32) | 1.59572 / 1.8 | — | 2.05088e-05 / 2.15733e-05 | no |
| curvature / Vort (8→16) | 1.69122 / 1.7 | — | 8.31519e-05 / 9.53425e-05 | no |
| curvature / Vort (16→32) | 1.04073 / 1.8 | — | 4.04185e-05 / 4.30907e-05 | no |
| curvature / Eion (8→16) | 1.66692 / 1.7 | — | 9.23384e-05 / 0.000100266 | no |
| curvature / Eion (16→32) | 1.71708 / 1.8 | — | 2.8086e-05 / 2.95534e-05 | no |

Measured runner time: 10442.3 s; 28 order checks fail. Recovery, error caps and diagnostics on/off checks pass.

CTest allows 3600 s for default and 18000 s for long (per-case limits
3600 / 10800 s). Grids, end time, error caps and order thresholds are retained.
Before/after timings include changed forcing-expression costs and concurrent
validation load; they are not isolated operator benchmarks.
