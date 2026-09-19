.. _sec-low-density-regions:

Low-density regions
===================

Very low density regions are numerically difficult in Hermes-3 because
several closures and transport terms require division by density or
temperature, while the physically relevant timescales also become very
short. This can lead to negative pressure, very large velocities,
non-finite source terms, or poor nonlinear solver robustness.

Hermes-3 therefore uses a collection of regularisations and robustness
options in low-density regions. These improve robustness, but they do
not remove the need for case-by-case tuning and they should not be
interpreted as a complete physical model of vacuum regions.


Current approach
----------------

The main ingredients are:

- Floors and soft floors. Density, temperature, and pressure can be
  limited using ``density_floor``, ``temperature_floor``, and derived
  pressure floors. In many places Hermes uses
  :math:`N_{lim} = \mathrm{softFloor}(N, N_{floor})` rather than a hard
  clamp, so the right-hand side remains differentiable.
- Low-density diffusion. Several components provide artificial
  parallel or perpendicular diffusion terms which become active where
  density or pressure is small, for example ``low_n_diffuse``,
  ``low_n_diffuse_perp``, ``low_T_diffuse_perp``, and
  ``low_p_diffuse_perp``.
- Equation-of-state regularisation. In some pressure-based models the
  evolved solver variable is treated as an internal-energy-like
  quantity based on :math:`N_{lim}`, while the physical pressure
  exposed to the rest of Hermes remains :math:`P = N T`.
- Optional freezing. The neutral fluid model ``neutral_mixed`` also
  provides ``freeze_low_density`` to suppress time evolution in very
  dilute regions while retaining sources and sinks.


Modified low-density equation of state
--------------------------------------

The most important recent change is that some pressure-based components
now distinguish between:

- the physical pressure, :math:`P = N T`
- the solver variable used to evolve internal energy in low-density cells

In these components the specific internal energy is regularised so that

.. math::

   e = C_V T \frac{N_{lim}}{N}

which means the evolved internal-energy-like quantity is proportional to

.. math::

   N e \propto N_{lim} T

rather than :math:`N T` in cells where the density floor is active.

This avoids some singular behaviour when :math:`N \rightarrow 0`, while
preserving the physical equation of state used by the rest of the model
away from the floor.


Affected components
-------------------

The following Hermes components currently use low-density regularisation
directly:

- ``evolve_density`` provides low-density diffusion coefficients that
  can be consumed by other components.
- ``evolve_pressure`` uses a limited-density internal energy in low-density
  cells and then maps back to physical pressure.
- ``evolve_momentum`` calculates velocity and several transport terms
  using a density floor to avoid singular velocities.
- ``neutral_mixed`` uses the same internal-energy idea for neutral
  pressure while also providing ``freeze_low_density``.
- ``neutral_full_velocity`` applies an analogous treatment in the full
  neutral velocity model.

The exact options available vary by component, so the component
documentation in :ref:`sec-equations` remains the authoritative source
for names and defaults.


Available options
-----------------

Commonly used low-density robustness options include:

``density_floor``
   Minimum density scale used when dividing by density or constructing
   limited density fields.

``temperature_floor``
   Temperature scale used by some low-temperature diffusion and sound-speed
   calculations.

``low_n_diffuse``
   Adds extra parallel diffusion in ``evolve_density`` to smooth very
   low density regions.

``low_n_diffuse_perp``
   Adds artificial perpendicular diffusion at low density.

``low_T_diffuse_perp``
   Adds artificial perpendicular diffusion at low temperature.

``low_p_diffuse_perp``
   Adds artificial perpendicular diffusion at low pressure.

``freeze_low_density``
   In ``neutral_mixed``, damps the time derivatives in very dilute
   regions while retaining explicit sources and sinks.

``damp_p_nt``
   In ``evolve_pressure``, damps the evolved pressure back toward
   :math:`N T` when pressure becomes negative or density falls below the
   floor.


Practical guidance
------------------

- Start by setting physically motivated source and boundary conditions.
  Low-density regularisation works best as a stabilisation layer, not as
  the primary mechanism that determines the solution.
- Increase floors and artificial diffusion gradually. Large values can
  substantially change the solution in detached or weakly populated
  regions.
- If a case fails only in dilute neutral regions, ``freeze_low_density``
  in ``neutral_mixed`` can be a useful diagnostic or steady-state aid.
- When changing low-density settings, inspect temperature, velocity, and
  source diagnostics in the same region rather than monitoring density
  alone.


Limitations and future work
---------------------------

Handling low-density regions remains an open problem in Hermes-3.
Current regularisations improve robustness, but they do not fully solve
all failure modes. In particular:

- floors modify the effective equations in dilute cells
- artificial diffusion can improve convergence at the cost of added
  dissipation
- different components still use different regularisation strategies
- some cases may require additional terms or new options beyond those
  currently implemented

Future extensions to low-density handling should be documented on this
page, with short component-specific summaries kept in
:ref:`sec-equations`.
