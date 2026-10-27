# Multispecies total-ion-energy balance

This test uses `evolve_vorticity_ion_energy` with a common ion temperature and
two ions with unequal masses and charges: `d+` has `(A,Z)=(2,1)` and `impurity`
has `(A,Z)=(3,2)`. Both evolve density and parallel momentum. Electrons use
`quasineutral`, `evolve_pressure` and `zero_current`; this algebraic current
closure adds no evolved electron kinetic reservoir.

Initial densities are `0.5*(1+0.1*sin(2*pi*x-z))` and
`0.25*(1-0.05*sin(2*pi*x-z))`, so composition and actual mass density vary
spatially. Initial momenta are 0.1 and -0.05, giving distinct velocities and
nonzero parallel energy. Initial aggregate ion energy is `1.2*Pe`;
pressure/potential recovery determines the common temperature after accounting
for kinetic energy.

The curvature slab has 16 by 16 physical X-Z cells and four periodic Y cells,
with Y-uniform fields and identity transformation. X/Z are periodic.
Its X-Z dimensions are 5 cm by 5 cm, `Rxy=1.5 m`, `Nnorm=2e18 m^-3`,
`Tnorm=5 eV` and `Bnorm=0.35 T`. Reference atomic mass is 1.75.
The input is [data/BOUT.inp](data/BOUT.inp); fixed criteria and source/boundary
provenance are in [criteria.json](criteria.json).

## Power balance and species weights

The independent ion and plasma inventories are

$$
E_i=\int\left[\frac32\sum_s n_sT_i
           +\sum_s\frac12 A_sn_sv_{\parallel s}^2+K_\perp\right]dV,\qquad
E=E_i+\frac32\int p_e\,dV.
$$

Saved `Eion` includes both ion kinetic reservoirs. The Python inventory compares
it with the partition above and adds electron thermal energy once. Its
finite-volume perpendicular operator uses

$$
\psi=\sum_s\frac{A_sn_sT_i}{\bar A Z_s},\qquad
\Phi=\phi+\psi,\qquad c=\frac{\bar A}{B^2}.
$$

The transported polarisation energy coefficient is reconstructed independently:

$$
E_{\rm adv}=\phi+\sum_s\frac{A_s}{Z_s\rho}
              \left(\frac52n_sT_i+K_{\parallel s}\right),\qquad
\rho=\sum_s A_sn_s.
$$

It uses actual mass density and the configured mass-density floor, while
`Kperp` uses the Boussinesq reference mass. The test requires unequal `A/Z`,
spatially varying composition/mass density, nonzero parallel energy and a
nonzero kinetic contribution to `E_adv`. Polarisation flux reconstruction and
direct species-current partition checks remain later work; coefficient coverage
alone does not validate those fluxes.

Only the prescribed evolution/closure components are active. Density and
pressure sources are explicitly zero; conduction, viscosity and extra diffusion
are disabled. Before supplying zero external and outward powers, the runner
validates these settings, periodic topology, saved X halos and Y uniformity.
It checks cumulative and interval plasma residuals,

$$
R(t)=E(t)-E(0)-\int_0^tQ_{\rm ext}\,dt'
                    +\int_0^tP_{\rm out}\,dt'.
$$

Nonzero curvature exchange closes the separate electron thermal budget through
snapshot reconstruction of `-Pe*Curlb_B.z*DDZ(phi)`, using saved normalisations
and unfiltered FFT DDZ. Enabled diagnostics cross-check this work, its integrated
electrical partner, ion kinetic/pressure diagnostics, zero external sources
and zero integrated energy-flux divergences.

## Cadence and acceptance

The default `ci` profile runs to normalised time 120 at output intervals 12,
6 and 3 with diagnostics enabled, plus an off run at interval 3.
`smoke` uses intervals 6 and 3 plus the off run. Pair histories must agree;
removing optional evolving outputs must preserve inventory and exchange analysis.

Saved times must be regular. The initial study gave electron residuals of
`1.80e-9`, `4.50e-10` and `1.12e-10` of the plasma reference energy,
decreasing about fourfold per interval halving. Plasma drift was below `4e-16`.
Parallel energy is about 0.25% of the plasma inventory. The exercised exchange
and flow signals exceed the fixed balance tolerances.

## Running and reports

From the repository root:

```sh
.venv/bin/python tests/integrated/2D-ion-energy-multispecies/runtest --exe build/hermes-3
.venv/bin/python tests/integrated/2D-ion-energy-multispecies/runtest --exe build/hermes-3 --profile smoke
env FI_SOCKETS_IFACE=en0 ctest --test-dir build -R '^2D-ion-energy-multispecies$' --output-on-failure
```

With the workspace environment activated, `./runtest` also works here.
`--output-root` selects an existing parent for isolated results and `--timeout`
limits each simulation. Inputs, hashes, logs, exit status and per-run
`energy-report.json` are retained. Reports contain powers in watts,
inventory/cumulative histories in joules, coverage and named assertions.
`result.json` contains this case's diagnostic-pair and cadence comparisons.

Tracked finding: XPPM with `bndry_flux=false` suppresses transport across the
radial periodic seam. It remains globally conservative in this case;
periodic-interface fidelity, parallel transport with Y gradients, general
sources and open-boundary fluxes need separate tests. See the
[analysis module](../../../tools/energy_analysis/README.md) and
[delivery plan](../ion-energy-analysis-plan.md).
