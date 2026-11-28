# Conservative polarisation transport

These options apply to `implementation = conservative` in
`[polarisation_drift]`, paired with `evolve_vorticity_ion_energy`.

```ini
[polarisation_drift]
implementation = conservative
polarisation_density_threshold = 0.02
polarisation_transport_weights = central
```

The density threshold is in units of `Nnorm`, like `density_floor`. Its
default is zero, disabling the limiter. For positive thresholds, the particle
flux is the polarisation current multiplied by
`s = max(rho, 0) / softFloor(rho, average_atomic_mass * threshold)`, where
`rho = sum(A_s N_s)` is the actual ion mass density before flooring.
The factor approaches zero at low density and one at high density. Each X-Z
face uses the current's upwind donor; physical walls always use the adjacent
boundary cell. The donor changes at zero current, keeping the flux continuous.

Density, parallel momentum, thermal energy and parallel kinetic energy follow
the limited particle flux. The current and its work flux `(phi + psi) Jpol`
retain their full value, preserving the perpendicular-energy exchange.
Particle charge transport can therefore differ from the polarisation current.
Optional poloidal transport applies the limiter to cell-centred particle
weights; its work weight remains unscaled.

`polarisation_transport_weights = central` preserves the arithmetic face
weights. `upwind` uses the same donor for composition and per-particle thermal
energy and momentum. Wall work always uses the local potential midpoint.

Four scalars are saved even with `diagnose = false`:

- `polarisation_charge_mismatch_l2`: volume-weighted RMS of
  `DivJpol + sum(Z_s polarisation_density_source_s)`, in A m⁻³.
- `polarisation_charge_mismatch_linf`: maximum absolute mismatch, in A m⁻³.
- `polarisation_charge_mismatch_rel`: RMS mismatch divided by RMS `DivJpol`;
  reported as zero when both norms vanish, or infinity when the mismatch is
  nonzero and `DivJpol` vanishes.
- `polarisation_limiter_volume_fraction`: volume fraction where `s < 0.5`.

With `diagnose = true`, `polarisation_charge_mismatch` saves the corresponding
3D field. Reductions use physical cells and run only at output times.
The wall energy and thermal weights include the density limiter, so multiplying
the saved wall energy weight by `polarisation_flux_x` reconstructs wall power.
The wall work weight retains the full local potential midpoint.

With `diagnose = true`, the signed energy-flux divergences are also saved as:

- `DivEpol_work`: the full-current work flux `(phi + psi) Jpol`.
- `DivEpol_transport`: thermal and parallel kinetic energy carried by the limited
  particle flux, using the selected central or upwind weights.
- `DivEpol_phi_work` and `DivEpol_psi_work`: the two contributions to work.

Each includes the physical-wall corrections and optional poloidal contribution.
Positive divergence removes ion energy. `DivEpol = DivEpol_work +
DivEpol_transport` and `DivEpol_work = DivEpol_phi_work + DivEpol_psi_work`,
up to roundoff. The psi work is already included in `DivEpol_work`.
These diagnostics retain the original arithmetic in the evolved energy equation.

For component source and momentum-work diagnostics, see
[ion-energy-source-diagnostics.md](ion-energy-source-diagnostics.md).

Energy and each species' particle wall losses are saved as
`balance_polarisation_drift_<quantity>_radial_polarisation_<region>`, where
`quantity` is `ion_energy` or `particles_<species>` and `region` is `core`,
`pf` or `outer`. The inner-X rank is core when its inner field line is
periodic in Y; otherwise it is PF. Positive values denote outward losses.
The previous `radial_polarisation_inner` energy entry is replaced by core
and PF entries. All three entries are present on every run, including zero
entries for absent walls.
