# Three-dimensional ion-energy slab MMS

One continuum manufactured solution extends the former `parallel_current_wall`
case. It has two ions, independent electron momentum, nonzero parallel current,
g13=1/4, g33=17/16, viscosity=0.02 and magnetic curvature=0.1. The metric is
positive definite with determinant one. The ion pressure potential retains
nonzero radial wall curvature and Neumann guards.

The total-potential fluctuation combines the first and third radial sine
harmonics in sheared slab coordinates. Its vorticity and normal derivative
vanish at both walls, matching the closed viscous-current faces. The pressure
potential retains nonzero wall curvature, while phi has Dirichlet wall values
and nonzero tangential flow.

Chi has Neumann AC and zero DC wall values. Its Y-mean Reynolds equation has
an analytic Fourier particular solution and homogeneous boundary correction.
The parallel current is derived by integrating
`dY(Jpar) = 2*laplacian(chi) - div(JReynolds) - div(Jdiamagnetic)` in Y.
The remaining RHS has first and second Y harmonics with periodic elementary
integrals.
The current includes the skew-metric and magnetic-curvature contributions.
Electron velocity uses analytic Dirichlet wall data for this derived current.
Viscosity remains outside chi and contributes vorticity diffusion and the
electron share of viscous heating. All forcing uses continuum derivatives;
numerical potential damping is enabled and is not cancelled by forcing.

The runner checks final-time evolved and recovered fields at t=0.02, starting
with 8 x 16 x 8 physical cells and doubling all dimensions. Spatial acceptance
uses the fixed central radial half, X/L in [0.25,0.75], with L2/L-infinity
order minima 1.4/0.9. Absolute field errors and global norms are diagnostic outputs. Recovery, current continuity,
source ownership, completed vorticity rate and nonzero current/Reynolds/
curvature/viscous signals are checked throughout the physical domain.
`--verify` adds tighter temporal integration and serial/two-rank agreement.
Committed inputs run without SymPy; regeneration and symbolic checks need it.

This study is not a routine CI test. CTest registration is explicitly opt-in:

```sh
cmake -S . -B build -DHERMES_ION_ENERGY_MMS=ON
ctest --test-dir build -R '^ion-energy-mms$' --output-on-failure
```

Alternatively run from the source tree:

```sh
python tests/integrated/ion-energy-mms/mms.py --check --check-files
python -m pytest tests/python/test_ion_energy_mms.py tests/python/test_reynolds_transport.py
python tests/integrated/ion-energy-mms/runtest --exe build/hermes-3 --levels 2 --verify --timeout 10800
python tests/integrated/ion-energy-mms/runtest --exe build/hermes-3 --plot --levels 3 --nproc 8 --timeout 43200
```

The opt-in CTest entry uses two spatial levels and temporal/MPI verification.
Additional refinement levels are available through `--levels`; larger 3D grids
benefit from more ranks and longer per-launch timeouts.

The default rank count is two. Use `--nproc`, `--mpiexec` and `--mpi-args` for
other launchers. Two levels are the default; plotting defaults to four.
`--output-root` retains results under a chosen parent directory. PDF/SVG figures,
CSV norms and `result.json` retain failures. Spatial failures stop dependent
verification; plot runs collect further spatial levels while retaining failure.
Execution and physical-contract failures always stop immediately. Time limits
are per launch; there are no separate grid-dependent output deadlines. Fine
coupled runs may take hours. Conservation and symbolic
identities do not establish passing coupled spatial convergence.
