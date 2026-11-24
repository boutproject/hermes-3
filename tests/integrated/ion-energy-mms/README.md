# Three-dimensional ion-energy slab MMS

One independent SymPy continuum solution has two ions, evolved electron
momentum, nonzero parallel current, g13=1/4, g33=17/16, viscosity=0.02 and
curvature=0.1. Sources use analytic derivatives, never discrete operators.
The metric is positive definite with determinant one. Potential has fixed
radial wall values; pressure potential has Neumann guards and nonzero wall
curvature. Vorticity and its normal derivative vanish at the walls, matching
closed viscous faces. Chi excludes viscosity and has Neumann AC/zero DC wall
data. Its analytic current satisfies the continuum elliptic identity.

`mms.py` generates inputs on demand from `base.inp`. `analyse.py` computes
volume-weighted errors against the continuum solution, using saved production
fields for kinetic energy and pressure recovery. The runner starts with
8x16x8 physical cells and doubles all dimensions. At t=0.02, orders must be
at least 1.4 in L2 and 0.9 in L-infinity over the fixed central radial half.
Global norms are also saved. `--verify` checks temporal errors below 5% of
spatial errors and serial/MPI agreement to absolute 1e-10. All runs check
recovery, the auxiliary solve residual and nonzero transport signals.

This is a long, optional study. Requires SymPy, boutdata, matplotlib and
xhermes on its `total-energy` branch:

```sh
python -m pytest tests/python/test_ion_energy_mms.py
python tests/integrated/ion-energy-mms/runtest --exe build/hermes-3 \
  --levels 2 --verify --output-root build/mms-study
```

Use a fresh output directory. `--levels` and `--nproc` control refinement and
parallelism. Norms and figures are saved as `errors.csv` and `convergence.pdf`.
For CTest, configure `-DHERMES_ION_ENERGY_MMS=ON`; its 10-hour timeout applies
to the entire study. Generated inputs, native outputs and launch logs are kept.
