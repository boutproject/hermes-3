# 2D energy test

Evolves vorticity, electron density, and electron and ion pressure in a 2D
(32 x 32) X-Z box with magnetic curvature. Includes kinematic viscosity and the
ion polarisation compression term. There is no parallel dynamics (`ny = 1`,
`MYG = 0`), no sources, and no sheath: apart from the radial potential boundary
relaxation, the box is energetically closed.

## Mechanism

The initial pressure perturbation, `1 + 0.1 * sin(2*pi*x - z)`, is unstable to
interchange modes. The curvature (`bxcvz = 1/R^2`) plays the role that gravity
plays in the Rayleigh-Taylor instability. The instability grows, drives
perpendicular E x B flows, and saturates. Kinematic viscosity then damps the
flows and returns their energy to the plasma as heat, so that the system relaxes
towards a stable state with a flattened pressure profile.

Unlike gravitational Rayleigh-Taylor, there is no separate potential energy
reservoir. The "potential energy" released by the instability is part of the
thermal energy itself. Curvature enters the model through two terms that are
exact energetic partners:

- the divergence of the diamagnetic current in the vorticity equation,
  $\nabla\cdot\mathbf{J}_d = \nabla\cdot\left[(p_e + p_i)\nabla\times(\mathbf{b}/B)\right]$;
- the compression of the E x B flow in the pressure equations,
  $-p_s\nabla\cdot\mathbf{v}_E = -p_s\nabla\phi\cdot\nabla\times(\mathbf{b}/B)$,
  applied by the `vorticity` component as $-\mathbf{J}_{d,s}\cdot\nabla\phi$.

The second term describes adiabatic expansion of flux tubes: in a torus, plasma
interchanged outward moves into larger flux-tube volume and cools. The constant
curvature in this slab represents that through $\nabla\cdot\mathbf{v}_E\neq 0$, even
though the slab metric is uniform. The instability therefore converts thermal
energy of both species into perpendicular kinetic energy; viscosity converts it
back into heat.

## Conserved energy

With the Oberbeck-Boussinesq approximation used here, the energy

$$
E = \int\left[\frac{3}{2}\left(p_e + p_i\right) + K_\perp\right]dV,
\qquad
K_\perp = \frac{1}{2}\frac{\bar{A}}{B^2}\left|\nabla_\perp\Phi\right|^2
$$

is conserved, where $\Phi = \phi + \hat{p}$ is the generalised potential,
$\hat{p} = \sum_s A_s p_s/(\bar{A}Z_s)$ is the pressure contribution to the
vorticity (`Pi_hat` in `src/vorticity.cxx`), and $\bar{A}$ is
`average_atomic_mass`. The energy exchange channels in this test are:

| Channel | Terms | Pairing |
|---|---|---|
| Curvature (interchange drive) | $\phi\nabla\cdot\mathbf{J}_d$ in vorticity; $-p_s\nabla\cdot\mathbf{v}_E$ in pressures | $\phi$ half of the diamagnetic current work |
| Polarisation compression | $\hat{p}\,\nabla\cdot\mathbf{J}_d$ in vorticity; $(A_s p_s/\bar{A}Z_s)\nabla\cdot\mathbf{J}_d$ in pressures (`polarisation_drift`) | $\hat{p}$ half of the diamagnetic current work |
| Reynolds stress | Nonlinear E x B advection of vorticity | Redistributes $K_\perp$; no net work in the continuum, *only for the full form* (see below) |
| Viscosity | $\nabla\cdot(\nu\nabla_\perp\omega)$ in vorticity; `kinematic_viscous_heating` $= -\nu\nabla_\perp\omega\cdot\nabla_\perp\Phi$ in pressures | Dissipation of $K_\perp$ into heat, distributed by mass density (almost all to ions) |
| Magnetic drift advection | $\nabla\cdot(\tfrac52 p_s\mathbf{v}_{mag,s})$, $\nabla\cdot(n\mathbf{v}_{mag,e})$ | Divergence; periodic in z, so no net contribution |
| Radial potential boundary | `phi_boundary_relax` | Does work through the boundary; not conserved (see below) |

In an exactly conserving model, therefore,

$$
\frac{3}{2}\int\left(p_e + p_i\right)dV\Big|_{\rm end}
= \frac{3}{2}\int\left(p_e + p_i\right)dV\Big|_{0} + K_\perp(0) - K_\perp(\rm end).
$$

The thermal energy first falls as the flows grow, then recovers as viscosity
heats the plasma. It ends at its initial value plus whatever perpendicular
kinetic energy has been removed overall. The initial $K_\perp$ is not zero: the
initial vorticity is zero, so $\Phi$ is harmonic with boundary values set by
$\hat{p}$, giving a combined E x B and diamagnetic flow at $t = 0$. It is small,
of order $(k\rho_s\,\delta p/p)^2\sim 10^{-5}$ of the thermal energy for these
parameters ($\delta p/p = 0.1$, 5 eV, 0.35 T, 5 cm box, $k\rho_s\approx 0.08$).

## Entropy

Energy is conserved, so the final state is not one of lower energy. What the
instability releases is the *available* energy: the difference between the
initial energy and the lowest energy reachable by ideal, entropy-preserving
rearrangements of the plasma, here flattening the pressure across the curvature
direction. Viscosity, and the numerical dissipation of the upwind E x B
advection, then convert this into heat irreversibly. The relaxed state is
approximately the state of maximum entropy at fixed total energy and particle
number. With entropy

$$
S = \sum_s\int n_s\ln\left(\frac{p_s^{3/2}}{n_s^{5/2}}\right)dV,
$$

the E x B advection and compression preserve the specific entropy
$p_s/n_s^{5/3}$ of each fluid element, so they rearrange entropy without
producing it. The viscous heating $Q_\nu$ produces entropy at the rate
$\int Q_\nu/T_i\,dV\geq 0$. The magnetic drift fluxes are reversible but do not
preserve specific entropy pointwise, so entropy is not a strict Lyapunov
function of the discrete model; it is monitored rather than asserted. Electrons
receive almost none of the viscous heating, so their entropy should change much
less than that of the ions.

## Expected results

- Electron thermal energy falls, because electrons lose energy to the flows by
  compression but receive almost none of the viscous heating.
- Ion thermal energy rises, because the ions receive the viscous heating.
- The total energy $E$ is constant, apart from the relaxation work at the radial
  boundaries and numerical dissipation.
- The total thermal energy changes by no more than $K_\perp(0) - K_\perp(\rm end)$.
  A larger rise indicates an energy source in the discrete model.
- Ion entropy increases.

## Known non-conservative terms

The following terms are present in this test as configured and are not expected
to conserve energy exactly. Measuring them is the purpose of the energy budget
diagnostics; the current runner checks only pressure-change signs.

1. **Simplified Reynolds stress.** `exb_advection_simplified` defaults to `true`,
   which advects $\omega$ with $\mathbf{v}_E$ only. With
   `diamagnetic_polarisation = true` the work done by this term is
   $\int\omega\,\mathbf{v}_E\cdot\nabla\hat{p}\,dV$, which is not a divergence. The full
   form (`exb_advection_simplified = false`) conserves $K_\perp$ in the continuum.
2. **Potential boundary relaxation.** `phi_boundary_relax = true` (default
   timescale $\tau = 10^{-4}$ s) relaxes the boundary potential towards zero radial
   gradient, doing work $P_{\rm relax}\simeq -(\Delta x/2\tau)\oint(\bar{A}/B^2)\,
   g^{xx}\partial_x\bar\phi\,\partial_x\bar\Phi\,dS$ through the radial boundaries.
3. **Viscous heating discretisation.** The viscous term in the vorticity equation
   is a finite-volume flux, while `kinematic_viscous_heating` is a cell-centred
   product of gradients; the two cancel only to the truncation error of the
   discretisation.
4. **Upwind E x B advection.** The limited (XPPM) advection of vorticity
   dissipates $K_\perp$ at small scales, and this energy is not returned as heat.
   The same limited advection of density and pressure is in flux form, so it
   conserves thermal energy.

## What the test checks

The `runtest` script currently checks the signs of the changes in electron,
ion and total pressure. The total-pressure check (total pressure higher at the
end than at the start) is consistent with energy conservation only if the rise
is no larger than $K_\perp(0) - K_\perp(\rm end)$, which is very small here.
Passing it does not distinguish physical heating from an energy leak.

It is planned to replace that check with quantitative criteria:

1. The relative change in total energy $E$, after subtracting the time-integrated
   boundary relaxation work, is below a stated tolerance.
2. The change in total thermal energy is bounded by $K_\perp(0) - K_\perp(\rm end)$
   plus that tolerance.
3. The electron thermal energy falls and the ion thermal energy rises (retained).
4. Ion entropy increases (monitored, not asserted).

The test is very low resolution, so it does not check convergence.

## Separate total-ion-energy tests

This test uses `vorticity`, `evolve_pressure` for both species, and
`polarisation_drift`. Tests of the aggregate `evolve_vorticity_ion_energy`
formulation have their own inputs, runners and power-balance explanations:

- [Uniform equilibrium](../2D-ion-energy-equilibrium/README.md).
- [Curvature exchange](../2D-ion-energy-curvature/README.md).
- [Unequal-mass/charge ion mixture](../2D-ion-energy-multispecies/README.md).

## Publication resolution study

The separate [publication runner](publication-study.md) compares matched
pressure and total-ion-energy runs across configurable spatial resolutions.
It produces state images, energy traces, convergence figures and numeric tables
for fixed or solver-managed relaxing potential boundaries.
