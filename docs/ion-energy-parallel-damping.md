# Numerical parallel potential damping

`evolve_vorticity_ion_energy:phi_dissipation` defaults to `false`. Enable it
explicitly in three-dimensional cases that need damping of parallel potential
variations, and include the `sound_speed` component in `[hermes] components`.
The generated `parallel_current` MMS enables it; regeneration leaves the
original baseline, skew, curvature and viscosity variants unchanged.

After recovery, scheduling and sheath boundary completion, finalisation uses
exactly the legacy vorticity operator and the communicated collective sound
speed:

```cpp
S_phi = -FV::Div_par(-phi, 0.0, sound_speed);
```

The default `fixflux=true` boundary treatment is retained: with zero advection
velocity the numerical flux at a physical Y target is zero. Internal and
periodic Y faces use MC reconstruction and the maximum collective sound speed
at the face. This damps the alternating-Y potential mode that collocated C2
parallel derivatives cannot resolve. A uniform-Y potential gives zero damping.
The electron `fastest_wave` is not substituted for `sound_speed`. An enabled
option without completed sound speed raises a configuration error.

The completed physical-cell vorticity rate is
`ion_current_rhs + DivJviscous + ion_phi_dissipation`. Encoded radial potential
guards and their fixed/free-slip/no-slip wall rates retain their original
meaning. `ion_phi_dissipation` is reset on every RHS evaluation, is zero in
nonphysical cells, and is saved when ion-energy diagnostics are enabled. Its
units/conversion match `ddt(Vort)` (`A m^-3`, `qe*Nnorm*Omega_ci`).
`ion_phi_dissipation_enabled` saves the resolved option. Component version 7
adds this optional term to the version-6 formulation; the restart layout is
unchanged. Analysis accepts versions 4 through 7 and reconstructs the damping
from potential, thermodynamic fields and geometry.

This term is a numerical vorticity rate, rather than a measured physical
parallel current or external power sink. It has no explicit ion/electron
heating or total-energy source. Evolved `Eion` retains its physical equation;
recovery repartitions the fixed total ion energy between thermal, parallel
kinetic and perpendicular kinetic energy. In the unforced periodic damping
probe, perpendicular kinetic loss is recovered as thermal energy. Monotonic
perpendicular kinetic loss is not an acceptance condition for arbitrary warm-ion
states with evolving walls and physical transport.

`ion_chi`, also exposed as the diagnostic alias `ddt(Phi)`, remains an auxiliary
approximation obtained from the selected current RHS. Viscosity and numerical
phi damping are excluded from that solve. Check full rate reconstruction after
subtracting both enabled terms; retain auxiliary continuity against the selected
current RHS without adding numerical damping.

The new parallel-current MMS treats chi absolute error as diagnostic-only. It
still saves absolute L2/Linfinity and interior errors, reference RMS and relative
RMS error, and requires the same spatial orders, temporal agreement, finite
errors and MPI agreement as before. Other field caps remain unchanged. MMS
forcing remains derived from the undamped continuum equations: it does not
cancel the grid-dependent numerical damping.

The equation and developer guides describe the current/energy treatment and
the Reynolds wall convergence limits. The focused periodic damping regression
is in `tests/integrated/ion-energy-parallel-damping/`.
