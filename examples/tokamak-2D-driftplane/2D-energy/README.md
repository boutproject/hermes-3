# Total-ion-energy resolution scan

Interchange relaxation in a 5 cm X-Z slab with curvature radius 1.5 m,
magnetic field 0.35 T, and density 2e18 m^-3. Ion total energy and electron
pressure evolve with conservative polarisation transport, viscosity and
classical diffusion (`custom_D = 0.1`). There are no prescribed sources or
sheaths. The radial transport faces are closed.

The pressure perturbation has constant buffers within 0.1 of either radial
wall and a quintic taper of width 0.1. These widths stay fixed across grids.
Native initialization sets `Eion = 1.5 Pe` and zero vorticity. This example
uses only the total-energy formulation.

From the repository root, with boutdata and xhermes (`total-energy` branch)
installed:

```sh
python examples/tokamak-2D-driftplane/2D-energy/run_scan.py build/energy-scan \
  --exe build/hermes-3 --levels 3 --base-nx 32 --base-nz 32
python examples/tokamak-2D-driftplane/2D-energy/analyse_scan.py build/energy-scan
```

Each level doubles X and Z. `--end-time`, `--outputs`, and `--nproc` control
launches; `--set SECTION:OPTION=VALUE` overrides input options. Use a fresh
output directory for each scan. The analysis runs independently and writes
energy histories, adjacent-grid electron thermal-energy differences, and particle
inventories to `convergence.pdf` and `convergence.csv`. Energy is `Eion +
1.5 Pe` integrated over physical cells in SI units; no flux or kinetic-energy
operator is reproduced in Python. For full power and particle balances,
see the nearby 3D `turbulence-ion-energy` example.
