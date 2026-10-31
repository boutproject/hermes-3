# DC/AC polarisation boundary test

This orthogonal 5 cm X-Z curvature slab has physical X faces, periodic Z,
singleton Y without guards, uniform hydrogen density and a smooth pressure
`p_e=p_i=1+0.1*cos(2*pi*z/Lz)` independent of X. There are no parallel flows,
sheaths, sources or viscosity. The runner builds native restart metadata,
then installs matched `phi=0`, FV `Vort=L(p_i)` and
`Eion=1.5*p_i+Kperp(p_i)` at time zero. This avoids an initial radial pressure
or potential boundary mismatch.

Initially the Reynolds current vanishes and the curvature forcing is a pure
AC mode independent of X. Its auxiliary potential must therefore also be
independent of X, with the exact FV Fourier amplitude. The test checks all
physical cells, including both boundary-adjacent rows. Before transport is
closed, the AC potential-flow normal flux must already be zero. A Dirichlet
AC inverse followed by removal of its boundary flux fails this check and
produces a density source proportional to inverse cell width.

Run using the repository environment:

```sh
.venv/bin/python tests/integrated/2D-ion-energy-polarisation-boundary/runtest \
  --exe build/hermes-3
```

The CI matrix refines X at fixed Z (16, 32, 64 by 16), with fixed and relaxing
phi boundaries. It also checks square refinement, a diagnostics-off pair,
a tighter-tolerance pilot and a two-rank radial decomposition. `--profile
smoke` uses the two coarsest X sizes and the diagnostics-off pair. `--sizes`
allows additional X resolutions, including 128, and `--end-time` sets the
normalised duration. Outputs remain regular. Each invocation writes isolated
inputs, matched initial states, logs, execution metadata and `result.json`.

The independent xhermes audit uses `bd.options` and the shared Python
finite-volume inversion with DC Dirichlet/AC Neumann conditions. It checks
partition and forward residuals, optional diagnostic divergences, closed
particle/energy faces, face Gauss identities and source-free plasma energy
and particle balance. Reports separate boundary/interior maxima and DC/AC
current corrections. Initial analytic checks provide the regression gate;
evolved vorticity amplitudes are reported without assuming they must vanish.

DC normal current remains permitted by the auxiliary Dirichlet condition and
is closed in particle and energy transport. This test isolates the AC source
in orthogonal geometry. Mixed-metric currents, DC forcing and species weights
are covered by the C++ and Python operator tests. The original publication
sine initialisation and its long-time nonlinear boundary behaviour are
separate checks.
