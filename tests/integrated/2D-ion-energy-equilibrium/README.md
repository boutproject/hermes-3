# Uniform total-ion-energy equilibrium

This test uses `evolve_vorticity_ion_energy` for aggregate ion energy and
`evolve_pressure` for electrons. Two ions have masses/charges `(2,1)` and
`(3,2)`, densities 0.5 and 0.25, and common normalised temperature 1.
Quasineutrality gives electron density 1. Initial electron pressure is 1 and
ion energy is 1.125. The slab has 16 by 16 physical X-Z cells, singleton Y,
and no curvature, gradients, sources or flows.

[data/BOUT.inp](data/BOUT.inp) defines the case, with `Nnorm=2e18 m^-3`,
`Tnorm=5 eV` and `Bnorm=0.35 T`. `nx=20` includes four radial guards.

## Energy balance

The saved ion energy must equal the independent common-temperature inventory,

$$
E_i = \int E_{\rm ion}\,dV
    = \int\left[\frac32\sum_s n_s T_i
                 + \sum_s\frac12 A_s n_s v_{\parallel s}^2
                 + K_\perp\right]dV.
$$

Both kinetic contributions vanish. Adding electron thermal energy gives the
plasma inventory. Uniform fields, no flow and no source components establish
zero external and outward powers; a radial potential boundary label alone
would not establish this for a nonuniform case. The runner verifies uniform
fields and zero kinetic energies before checking constant energy.

Each actual diagnostic on/off run has four output intervals of 10 normalised
time. All records must pass ion partition, forward vorticity, species discovery
and optional kinetic/pressure diagnostic checks. State and inventory histories
must agree between runs. Removing all optional time-dependent fields from the
enabled dataset must preserve the inventory. Fixed limits in `runtest` require
relative energy drift and normalised partition/diagnostic errors below `1e-10`.

## Running and reports

From the repository root, use the workspace environment:

```sh
.venv/bin/python tests/integrated/2D-ion-energy-equilibrium/runtest --exe build/hermes-3
env FI_SOCKETS_IFACE=en0 ctest --test-dir build -R '^2D-ion-energy-equilibrium$' --output-on-failure
```

With the environment activated, `./runtest` also works from this directory.
`--output-root` selects an existing parent for a new isolated results
directory; `--timeout` limits each simulation. Inputs, run logs,
`execution.json`, `energy-report.json` and `result.json` are retained.

See the [analysis module](../../../tools/energy_analysis/README.md) for the
data contract and operators, and the
[delivery plan](../ion-energy-analysis-plan.md) for further coverage.
