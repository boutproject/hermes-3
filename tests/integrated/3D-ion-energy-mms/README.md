# Ion total-energy MMS in a 3D slab

This test uses the method of manufactured solutions (MMS) to check the
coupled two-ion model with a common ion temperature. It evolves `Ni`, `Nj`,
`NVi`, `NVj`, `Pe`, `Vort` and `Eion`, and independently checks the recovered
total and species ion pressures, `Tion`, `phi` and `psi`. Optional diagnostics
also allow the auxiliary `chi` to be compared with its analytic reference.

The normalised slab has `Lx=Ly=Lz=10`, identity metrics, `B=J=1`, zero curvature
and no sheaths. X has physical boundaries: homogeneous Dirichlet `phi` and
Neumann recovered pressure. Y and Z are periodic. Ion mass/charge pairs are
`(2,1)` and `(3,2)`, with quasineutral electrons and zero total parallel
current. Both ions retain density and momentum evolution; electrons retain
pressure evolution. Collisions, conduction and extra diffusion are absent.

`manufactured.py` differentiates smooth, strictly positive primitive profiles
in physical normalised coordinates. It derives vorticity and energy from the
continuum potential and kinetic partition. The reference is independent of
mesh spacing, including at time zero: recovery then has its own spatial
truncation error. MMS forcing is `dq_exact/dt - F_physical`, added once by
BOUT after the model RHS. It is not a physical species source and receives
no additional kinetic-energy conversion.

The auxiliary potential solves `L(chi)=F_Vort` before the manufactured source
is added, rather than `L(chi)=dq_exact/dt`. Its analytic cosine coefficients
satisfy both zero Dirichlet values and zero normal derivatives at X faces.
This satisfies the auxiliary midpoint Dirichlet DC / radial Neumann AC
current solve and closed normal polarisation transport. The generator checks these
identities symbolically and reports the magnitudes of individual RHS terms.

Run from the repository root using the workspace environment:

```sh
.venv/bin/python tests/integrated/3D-ion-energy-mms/mms.py --check
.venv/bin/python tests/integrated/3D-ion-energy-mms/runtest \
  --exe build/hermes-3 --output-root build/mms-scan
```

The default and registered scan is `8^3, 16^3, 32^3`, to normalised time 0.1,
with one output interval. This keeps regression testing bounded while measuring
the continuum recovery error at time zero and the evolved error later.
Use `--sizes`, `--end-time` and `--outputs` for longer, finer studies. Larger
studies need an appropriate `--timeout` (seconds per case): an exploratory
64^3 baseline case exceeded the default 900-second limit, so refinement beyond
32^3 has not been validated. A short exploratory run can use `--sizes 8 16 --end-time 0.1 --outputs 1
--no-criteria`. Results remain in isolated case directories, with inputs,
resolved settings, logs and execution metadata. A nonzero launcher return
code is accepted only if every rank has written all expected fields and
times; analysis checks the output axis again.

Additional controls are:

- `--diagnose`: save and check the auxiliary potential.
- `--tolerance-check`: repeat the finest grid with tighter time/recovery
  tolerances and half the maximum timestep; field differences must be below
  10% of their spatial L2 errors.
- `--mpi-checks`: repeat the medium grid with X, Y and 2x2 X-Y decompositions.
  The main scan is serial; `--mpiexec` and `--mpi-args` configure the launcher.
- `--axis x`, `y` or `z`: refine only that direction, holding others at the
  first requested size. This diagnoses anisotropy; use `--no-criteria` when
  errors in the fixed directions dominate.
- `--nproc`/`--nxpe`: run the main scan in a chosen MPI decomposition.
- `--resume PATH`: reuse completed outputs with matching inputs/settings.
- `--analyse-only PATH`: regenerate reports without starting simulations.
- `--expression-sources`: use the ordinary BOUT expression evaluator; this
  also works in builds with `HERMES_TESTS=OFF`.

Test builds normally use a compiled evaluator generated from exactly the same
SymPy sources as `data/BOUT.inp`. Before overriding BOUT's source generators,
it compares the compiled and BOUT expression values at four interior/boundary
and time-dependent points on every rank. A shared cache evaluates all seven
sources once per cell and RHS time. BOUT still adds the MMS forcing and computes
native errors. The original analytic source expressions are preserved in
resolved settings under `test_ion_energy_mms`. Physical vorticity forcing is
masked outside the radial interior in both source backends: BOUT adds MMS
sources after the component RHS, and radial `Vort` guards now store fixed phi
midpoints. The audit checks Neumann physical vorticity guards and the
mean-Neumann, zero-AC-midpoint pressure/potential guards from Plan C.

CVODE uses the Adams method with fixed-point iteration. The BOUT change in this
implementation prevents calls to linear-solver setup/statistics routines in
that mode, and its solver regression covers the setup. Default time tolerances
are `atol=1e-10`, `rtol=1e-8`, with maximum timestep 0.05. Potential recovery
uses `phi_atol=1`, `phi_rtol=1e-13`, `vort_atol=1e-12`: `phi_atol` is the scale
in the denominator of the relative update test, so its product with `phi_rtol`
sets the absolute update tolerance near potential zeros.

The analysis opens an xhermes dataset with `unnormalise=False` and uses
`bd.options` for resolved settings. The raw multiline git diff in `[run]`
is omitted from an analysis copy of the settings because it can confuse the
options reader; all resolved physics options are retained. It constructs the physical cell-centre
coordinates independently and cross-checks native `E_<variable>` errors
against direct Python evaluation. A mismatch fails even exploratory runs.
Volume-weighted RMS L2, Linf, relative errors, adjacent-grid orders, error
locations, boundary face errors and boundary-adjacent cell errors are saved
for every output time. Guards and duplicate processor cells are excluded;
all physical cells, including extrema and cells beside X faces, are retained.
Dirichlet phi fixes the gauge: no mean subtraction is performed.

`errors.csv` and `errors.json` contain the scan tables. `convergence.pdf` and
`convergence.png` plot final errors. Each case also reports finite-volume
vorticity and energy-partition residuals, reconstructed with the shared Python
operators. Those residuals measure recovery convergence, whereas comparison
with continuum `phi` and pressure measures spatial accuracy. MMS explicitly
forces the system, so global energy conservation is not an acceptance gate.

Production MC limiters in X-Z and the configured parallel limiter remain
active. Smooth solutions still have extrema; Linf order need not be two.
The final-time orders measured between 16^3 and 32^3 are:

| Field | Baseline L2 / Linf | Skew L2 / Linf | Viscosity L2 / Linf |
| --- | --- | --- | --- |
| Ni | 1.73 / 1.22 | 1.69 / 1.33 | 1.69 / 1.29 |
| Nj | 1.78 / 1.16 | 1.75 / 1.23 | 1.75 / 1.27 |
| NVi | 1.81 / 1.36 | 1.81 / 1.40 | 1.81 / 1.40 |
| NVj | 1.85 / 1.23 | 1.85 / 1.32 | 1.85 / 1.31 |
| Pe | 1.67 / 1.14 | 1.67 / 1.20 | 1.67 / 1.19 |
| Vort | 1.50 / 1.11 | 2.13 / 2.18 | 1.91 / 1.66 |
| Eion | 1.91 / 1.70 | 1.83 / 1.30 | 1.83 / 1.30 |
| Pi | 1.96 / 1.74 | 1.91 / 1.62 | 1.91 / 1.64 |
| phi | 2.00 / 1.92 | 2.02 / 1.89 | 2.01 / 1.89 |

`criteria.json` gates both norms on the finest pair, with explicit recovered-field
thresholds and margin below the measured orders. It requires three resolutions.
In the baseline tolerance audit, the largest change in L2 was below 0.006% of
the spatial error; serial/X/Y/2x2 decompositions agreed within 4e-14 in Linf.
The corresponding tolerance changes for skew/viscosity were below 0.13% of
the spatial errors, and their MPI comparisons agreed within 3e-15.

The baseline's derived `DelpPhi_2B2` stencil retains the model's Neumann guards,
although the continuum normal derivative is nonzero. This does not dominate
the measured full-domain vorticity error: at 32^3 and t=0.1, the boundary-cell
Linf is 2.13e-5 versus 2.14e-4 over the domain, whose maximum lies at an interior
extremum. Boundary-cell L2/Linf and face-condition errors are retained explicitly
in the reports. The bounded scan does not establish the asymptotic behaviour of
this stencil at higher resolution.

The [skew-metric](../3D-ion-energy-mms-skew/README.md) and
[viscosity](../3D-ion-energy-mms-viscosity/README.md) tests have separate inputs,
READMEs and runners. Generate their compatible profiles with `mms.py --all
--check`; changing metric or viscosity options in the baseline alone does not
preserve its manufactured boundary problem.
# Auxiliary inversion accuracy

`--corrections N` selects the fixed FV correction count after the initial
spectral chi solve (default two) and records it in the scan manifest. Repeat
with `--corrections 6` to measure sensitivity without changing continuum MMS
sources or pressure/potential recovery tolerances. With `--diagnose`, chi's
continuum errors are convergence measurements; a nonzero discrete auxiliary
residual is distinct from pressure/potential recovery and energy partition
errors.
