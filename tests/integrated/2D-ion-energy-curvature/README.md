# Curvature-driven total-ion-energy balance

This test evolves electron density and pressure, one quasineutral ion with
mass/charge `(1,1)`, and aggregate ion energy with
`evolve_vorticity_ion_energy`. The pressure perturbation
`1+0.1*sin(2*pi*x-z)` drives interchange motion through constant magnetic
curvature `bxcvz=1/Rxy^2`. Initial ion energy is `1.5*Pe`. The test measures
early electron/ion energy transfer, before a long saturated interchange run.

The slab has 16 by 16 physical X-Z cells, singleton Y, periodic X/Z and closed
core Y connections. Its X-Z dimensions are 5 cm by 5 cm, with `Rxy=1.5 m`,
`Nnorm=2e18 m^-3`, `Tnorm=5 eV` and `Bnorm=0.35 T`. There is no parallel
dynamics, conduction, viscosity or prescribed heating. The input is
[data/BOUT.inp](data/BOUT.inp); fixed acceptance limits and source/boundary
provenance are in [criteria.json](criteria.json).

## Power balance

The plasma inventory is

$$
E=\int\left(E_{\rm ion}+\frac32p_e\right)dV,\qquad
E_{\rm ion}=\frac32n_iT_i+K_\perp.
$$

The independent `Kperp` calculation uses the finite-volume face/corner operator
on `Phi=phi+psi`, where `psi=p_i` for this ion and reference mass. Reconstructed
ion energy and forward vorticity must agree with saved fields at every record,
including startup.

The controlled case has zero external and outward powers. Before supplying
them, the runner validates the component whitelist, explicit zero sources,
saved periodic topology and X halos, Y uniformity, and disabled dissipation.
It checks cumulative and interval residuals,

$$
R(t)=E(t)-E(0)-\int_0^tQ_{\rm ext}\,dt'
                    +\int_0^tP_{\rm out}\,dt'.
$$

The electron thermal budget exercises a nonzero internal exchange:

$$
W_e=-p_e\,{\rm Curlb_B}^{z}\,\partial_z\phi,\qquad
Q_e=\int W_e\,dV,\qquad
\Delta E_e=\int_0^tQ_e\,dt'.
$$

`electron_curvature_work` reconstructs this from ordinary snapshots, scalar
curvature, saved normalisations and unfiltered FFT DDZ. This spectral derivative
is distinct from the finite-volume perpendicular-energy discretisation.
Analysis works with diagnostics disabled. Enabled runs compare
`electron_exb_work` pointwise, `phi_Div_Je` after integration, and ion
kinetic/pressure diagnostics. They also check zero external source diagnostics
and zero integrals of each `DivE*` flux divergence.

## Cadence and acceptance

All runs end at normalised time 120. The default `ci` profile uses output
intervals 12, 6 and 3 with diagnostics enabled, plus an actual off run at interval
3. `smoke` uses intervals 6 and 3 plus the off run. Pair histories must agree;
stripping optional time-dependent outputs must preserve inventory and exchange
analysis.

Saved times must be regular. Trapezoidal exchange residuals must decrease with
cadence refinement; the third cadence tests the smooth-power assumption.
The initial study found residuals of `3.34e-9`, `8.37e-10` and `2.09e-10`
of the plasma reference energy, about fourfold reduction per interval halving.
Plasma drift was below `4e-16`. Nonzero exchange and perpendicular energy must
exceed the fixed balance tolerances comfortably.

## Running and reports

From the repository root:

```sh
.venv/bin/python tests/integrated/2D-ion-energy-curvature/runtest --exe build/hermes-3
.venv/bin/python tests/integrated/2D-ion-energy-curvature/runtest --exe build/hermes-3 --profile smoke
env FI_SOCKETS_IFACE=en0 ctest --test-dir build -R '^2D-ion-energy-curvature$' --output-on-failure
```

With the workspace environment activated, `./runtest` also works here.
`--output-root` and `--timeout` select the results parent and per-run limit.
Each run retains inputs, hashes, logs, exit status and `energy-report.json`
with powers in watts, inventories/cumulative contributions in joules, coverage
and assertions. `result.json` holds this case's pair and cadence comparisons.
`analysis.py` provides the generic saved-run analysis CLI.

Tracked finding: XPPM with `bndry_flux=false` suppresses transport across the
radial periodic seam through `firstX/lastX` checks. This case remains globally
conservative; correct periodic-interface transport is a separate check.

[reference/BOUT.inp](reference/BOUT.inp) preserves the older `data-Eion` input.
It has radial boundaries and is a standalone long-running reference, rather
than a closed-domain conservation gate. `runtest` uses only `data/BOUT.inp`.
See the [analysis module](../../../tools/energy_analysis/README.md) and
[delivery plan](../ion-energy-analysis-plan.md) for further flux/boundary work.
