# COMPASS potential relaxation and sheath fluxes

This axisymmetric case uses the existing 36 × 48 COMPASS grid, one toroidal
cell and six MPI ranks. It contains a closed core, scrape-off layer, private
flux region and two sheath targets. Two ions with unequal mass-to-charge
ratios share one temperature. Electron pressure uses `evolve_pressure`.
Imposed particle and heating sources are explicitly zero; parallel flow and
sheath losses remain enabled. The radial potential midpoints relax on a
1 microsecond timescale, with a physical-cell Y average on the connected core.

Run using the workspace environment:

```sh
.venv/bin/python tests/integrated/tokamak-ion-energy-relax-boundary/runtest --exe build/hermes-3
```

The runner stages the grid from `examples/tokamak-2D/` in each isolated run.
Build-tree fixtures receive a copy through CMake. `MPIEXEC` can select a
launcher; the default is `mpiexec -n 6`. The default MPICH launcher supports
six ranks on this small grid without an additional oversubscription flag.

The independent analysis reads each rank with its processor halos. It checks
the ion energy partition, X-Z forward vorticity operator, midpoint relaxation
ODE, sheath thermal source assembly and parallel ion enthalpy/kinetic boundary
power. Parallel face fluxes must cancel across ordinary processor interfaces,
the closed-core connection and the private-flux connection. The sheath
reconstruction uses completed boundary values reconstructed from interiors,
since saved solver guards need not contain the final sheath state.

With diagnostics enabled, the saved auxiliary `ddt(Phi)` supplies χ to an
independent X-Z face reconstruction. It must reproduce `DivJpol` and `DivEpol`
with exactly zero radial particle, momentum, enthalpy, parallel kinetic and
polarisation electrical faces, while its unmasked normal current is nonzero.
Axisymmetry and disabled poloidal flows make these polarisation Y fluxes zero
at both targets. This geometry exercises topology and boundary powers; the
slab and unit tests exercise nontrivial Z structure and mixed-metric fluxes.

The CI matrix runs three regular output cadences and a diagnostics-off pair.
`--profile smoke` uses the two finest cadences. Both modes reconstruct
inventories, relaxation and sheath/parallel powers from ordinary outputs.
Diagnostics provide additional source/flux cross-checks and a sampled ion RHS
balance. Reports explicitly mark general tokamak power coverage incomplete:
electrical/diamagnetic boundary work, the independent auxiliary-current RHS
and the full electron power/kinetic exchange audit are not yet reconstructed.
The diagnostic RHS balance does not establish a complete independent power
audit. No absent contribution is treated as zero.

The audit preserves the current simple-sheath geometry and signed particle
diagnostic conventions. In particular, its lower ion source uses the next
interior metric when calculating sheath area; the parallel finite-volume
operator averages cell-face areas. The reports distinguish thermal sheath
losses from ion advected boundary power and do not compensate for their
different conventions with an invented source.
