.. _sec-boundary_conditions:

Boundary conditions
-------------------

Boundary conditions in Hermes-3 are applied at several levels:

Fundamental boundary conditions
   These are implemented in BOUT++ and typically default to dirichlet(0)
   with some exceptions.

Component-set boundary conditions
  Individual components, such as ``sheath_boundary``, will override
  the fundamental BCs, for example to set a Bohm condition and 
  extrapolate to the wall.

Generally speaking, it shouldn't be necessary to change any fundamental
boundary conditions, as they are handled by the components. A notable
exception to this is at the radial boundaries, which default to ``Neumann``
(zero gradient), preventing any cross-field transport. This can be
set to ``dirichlet(value)`` to set a constant value for a particular quantity.
The outer radial boundary can be set to ``decaylength(value)`` to allow 
losses according to a far-SOL decay length, see below for details.
Eventually, a component will be added to handle this at a higher level.

Fundamental boundary conditions
~~~~~~~~~~~~~~~

BOUT++ provides a number of fundamental boundary conditions including:

- dirichlet(x): boundary set to constant value of `x`
- neumann: boundary set to zero gradient
- free_o2: boundary set by linear extrapolation (using 2 points)
- free_o3: boundary set by quadratic extrapolation (using 3 points)

These can be set on different parts of the domain using the keywords
`core`, `sol`, `pf`, `lower_target`, `upper_target`, `xin`, `xout`, `yup`, `ydown` and `bndry_all`.

The boundary conditions can also be applied over a finite width as well as relaxed over a specified timescale.

These boundary conditions are implemented in BOUT++, and therefore have no access to
the normalisations within Hermes-3 and so must be used in normalised units.
Please see the `BOUT++ documentation
<https://bout-dev.readthedocs.io/en/latest/user_docs/boundary_options.html>`_ for more detail, 
including the full list of boundary conditions and more guidance on their use.
In case the documentation is incomplete or insufficient, please refer to the 
`BOUT++ boundary condition code
<https://github.com/boutproject/BOUT-dev/blob/cbd197e78f7d52721188badfd7c38a0a540a82bd/src/mesh/boundary_standard.cxx>`_
.


Currently, there is only one additional simple boundary condition implemented in Hermes-3.
`decaylength(x)` sets the boundary according to a user-set radial decay length. 
This is a commonly used setting for plasma density and pressure in the tokamak SOL boundary in 2D and 3D but is not applicable in 1D.
Note that this must be provided in normalised units just like the BOUT++ simple boundary conditions.

The below example for a 2D tokamak simulation sets the electron density to a constant value of :math:`1e20^{-3}` in the core and
sets a decay length of 3mm in the SOL and PFR regions, while setting the remaining boundaries to `neumann`.
Example settings of the fundamental normalisation factors and the calculation of the derived ones is provided
in the `hermes` component which can be accessed by using the `hermes:` prefix in any other component in the input file.

.. code-block:: ini

   [hermes]
   Nnorm = 1e17  # Reference density [m^-3]
   Bnorm = 1   # Reference magnetic field [T]
   Tnorm = 100   # Reference temperature [eV]
   qe = 1.60218e-19   # Electron charge
   Mp = 1.67262e-27   # Proton mass
   Cs0 = sqrt(qe * Tnorm / Mp)   # Reference speed [m/s]
   Omega_ci = qe * Bnorm / Mp   # Reference frequency [1/s]
   rho_s0 = Cs0 / Omega_ci   # Refence length [m]

   [Ne]
   bndry_core = dirichlet(1e20 / hermes:Nnorm)
   bndry_sol = decaylength(0.003 / hermes:rho_s0)
   bndry_pf = decaylength(0.003 / hermes:rho_s0)
   bndry_all = neumann()


Sheath
~~~~~~~~~~~~~~~
Hermes-3 includes additional boundary conditions whose complexity requires their implementation
as components. They may overwrite simple boundary conditions and must be set in the same way as other components.

.. _sheath_boundary_simple:

sheath_boundary_simple
^^^^^^^^^^^^^^^

This is a top-level component which determines the conditions and sources at the divertor target. 
First, density, temperature and pressure are extrapolated into the target boundary.
The extrapolation method for each can be user set, e.g. `density_boundary_mode`. At the moment, 
the available modes are:

- 0: LimitFree
   An exponential extrapolation for decreasing quantities and a Neumann boundary for increasing
   quantities. It is inconsistent between increasing and decreasing values and is not recommended
   unless you have a reason to use it - it's legacy behaviour.

- 1: ExponentialFree
   An exponential extrapolation. It is more consistent than LimitFree and has the advantage of 
   inherently preventing negative values at the target. This is the default and is recommended for most cases.
   It is defined as :math:`guard = (last)^2 / previous`

- 2: LinearFree
   A linear extrapolation. It can lead to negative values at the target. However, it is the most 
   consistent (the linear extrapolation is what second order differencing reduces to at the wall) and has been shown to 
   reduce "zigzags" or "squiggles" near the target which are common in cell centered codes. Use only
   if you have a particular reason to care about this.
   It's defined as :math:`guard = 2 * last - previous`.

In the above definitions, `last`, `previous` and `guard` refer to the final domain cell, the penultimate 
domain cell and the guard cell respectively. The value at the target is defined to be 
an interpolation between the last and guard cells, i.e:

.. math::
   \begin{aligned}
   target = (last + guard)/2
   \end{aligned}

After the initial extrapolation, the sheath velocity is set to greater or equal to Bohm speed as according
to Stangeby, eq. 2.55b, with temperature in eV:

.. math::
   \begin{aligned}
   v_{i}^{sheath} \geq [(e T_{e} + \gamma e T_{i})/m_{i}]^{1/2}
   \end{aligned}

where :math:`\gamma` is the ion polytropic coefficient, which is set to 1 by default as per SOLPS-ITER.
The electron velocity is calculated from the potential:

.. math::
   \begin{aligned}
   \phi^{sheath} &= T_{e} \  ln \biggl[ \sqrt{ T_{e} / (2 \pi m_e) \cdot (1 - G_{e}) \cdot n_{e} / \Gamma_{i}^{tot}} \biggr]   \\ 
   v_{e}^{sheath} &= -\sqrt{ T_{e} / 2 \pi m_e} \cdot (1-G_{e}) \cdot e^{(- \frac{\phi_{sheath} - \phi_{wall}} {T_{e}})}
   \end{aligned}

Where :math:`\Gamma_{i}^{tot}` is the total ion particle flux across all species, :math:`G_{e}` is the electron secondary
emission coefficient (default 1) and :math:`\phi_{wall}` is the wall potential (default 0). Both can be user-set through the
options `secondary_electron_coef` and `wall_potential`, respectively.


Hermes-3 allows the pressure and momentum equations to advect internal and 
kinetic energy out of the sheath, but disables the conduction term, so before applying any sheath boundary condition,
the "underlying" sheath heat flux is:

.. math::
   \begin{aligned}
   q_{i}^{sheath} &= n_{i} v_{i} (\frac{5}{2} e T_{i}n_{i} + \frac{1}{2} m_i v_{i}^{2})

   q_{e}^{sheath} &= n_{e} v_{e} (\frac{5}{2} e T_{e}n_{e} + \frac{1}{2} m_e v_{e}^{2})
   \end{aligned}

With both of the above definitions following Stangeby (eqns. 9.61 and 9.63) with the addition of electron kinetic
energy, and where each variable is evaluated at the sheath (target). Continuing from Stangeby (eqs. 2.89 and 2.94), 
the above can be represented as the internal energy multiplied by the sheath heat transfer coefficient:

.. math::
   \begin{aligned}
   q_{i}^{sheath} &= n_{i} v_{i} \gamma_{i} (e T_{i}n_{i})
   q_{e}^{sheath} &= n_{e} v_{e} \gamma_{e} (e T_{e}n_{e})
   \end{aligned}

Where :math:`\gamma_{i}` and :math:`\gamma_{e}` are the **total** ion and electron sheath heat transfer coefficients 
and are user set through the options `gamma_i` and `gamma_e`. The easiest way to calculate the sheath heat flux leaving the model
is to use this form of the equations.

Assuming :math:`T_e = T_i`, the "underlying" heat flux corresponds to :math:`\gamma_{i} = 3.5` and :math:`\gamma_{e} = 3.5`.
In order to facilitate a user-set total gamma, `sheath_boundary_simple` creates a heat sink (or source) 
based on the difference between the required and "underlying" coefficient:

.. math::
   \begin{aligned}
   q_{i}^{additional} &= \gamma_i T_i n_i v_i - (2.5 T_{i} + \frac{1}{2} m_i v_{i}^{2}) n_{i} v_{i}

   q_{e}^{additional} &= \gamma_e T_e n_e v_e - (2.5 T_{e} + \frac{1}{2} m_e v_{e}^{2}) n_{e} v_{e}
   \end{aligned}

The sheath ion particle flux is facilitated through the underlying density equation advecting density out of the domain:

.. math::
   \begin{aligned}
   \Gamma_{i}^{sheath}&= n_{i} v_{i}
   \end{aligned}

If recycling is enabled, a corresponding recycled neutral particle source is set up in the `recycling` component.

**Usage and other options**

To enable the boundary condition, add `sheath_boundary_simple` to the list of components. The settings are
accessed through the component header. Use the `lower_y` and `upper_y` flags to enable the boundary on each 
end of the domain, where `lower` and `upper` refer to the start and end of the poloidal index space, respectively:

.. code-block:: ini

   [hermes]
   components = d+, sheath_boundary_simple

   [d+]
   type = noflow_boundary

   noflow_lower_y = true   # This is the default
   noflow_upper_y = false  # Turn off no-flow at upper y for d+ species

   [sheath_boundary_simple]
   lower_y = false         # Turn off sheath lower boundary for all species
   upper_y = true

It can be useful to run the code without neutrals/recycling in order to simplify the physics, e.g. for debugging.
However, just disabling `recycling` would result in mach flows throughout the whole domain due to the lack of the
neutral source. To avoid this, you can set `no_flow = true` under `sheath_boundary_simple`. This will set the ion 
velocity to zero for the particle flux but will keep it at the :math: `v_i \geq c_{bohm}` condition for the heat flux.

By default, the Bohm condition is imposed on the target by the Lax flux. This allows the code to have a small amount 
of slack, resulting in a not-perfectly-exact setting but a smoother and more stable solution. For debugging, you can
disable this behaviour and fix the Bohm condition explicitly. This can be done by setting `fix_momentum_boundary_flux` 
to `true` in the `evolve_pressure` component. Note that this has been observed to increase numerical oscillations near
the boundary and is not recommended.

.. _sheath_boundary:

sheath_boundary
^^^^^^^^^^^^^^^

This component is required to calculate correct sheath heat transfer coefficients considering multiple main ions
based on Tskhakaya 2005. As this component is more complex, the development may lag behind `sheath_boundary_simple`.

.. _sheath_boundary_insulating:

sheath_boundary_insulating
^^^^^^^^^^^^^^^

WIP

.. _noflow_boundary:

noflow_boundary
^^^^^^^^^^^^^^^

WIP

Recycling
~~~~~~~~~

This component calculates the flux of a species into a boundary
due to recycling of flow out of the boundary of another species.

The boundary fluxes might be set by sheath boundary conditions,
which potentially depend on the density and temperature of all species.
Recycling therefore can't be calculated until all species boundary conditions
have been set. It is therefore expected that this component is a top-level
component (i.e. in the `Hermes` section) which comes after boundary conditions are set.

Recycling has been implemented at the target, the SOL edge and the PFR edge.
Each is off by default and must be activated with a separate flag. Each can be 
assigned a separate recycle multiplier and recycle energy. 

Configuring thermal recycling
^^^^^^^^^^^^^^^

A simple and commonly used way to model recycling is to assume it is fully thermal,
i.e. that every incident ion recombines into a neutral molecule and thermalises with the surface 
before becoming re-emitted. Hermes-3 does not yet have a hydrogenic molecule model, and so 
the molecules are assumed to instantly dissociate at the Franck-Condon dissociation temperature of 3.5eV.

In order to set this up, the chosen species must feature an outflow through the boundary - any cells
with an inflow have their recycling source set to zero. If a sheath boundary condition
is enabled, then this is automatically satisfied at the target through the Bohm condition.
If it is not enabled, then the target boundary must be set to `free_o2`, `free_o3` or `decaylength` to 
allow an outflow. 

The recycling component has a `species` option, that is a list of species
to recycle. For each of the species in that list, `recycling` will look in
the corresponding section for the options `recycle_as`, `recycle_multiplier`
and `recycle_energy` for each of the three implemented boundaries. Note that 
the resulting recycling source is a simple
multiplication of the outgoing species flow and the multiplier factor.
This means that recycling `d+` ions into `d2` molecules would require a multiplier 
of 0.5 to maintain a particle balance in the simulation.

For example, recycling `d+` ions into `d` atoms with a recycling fraction
of 0.95 at the target and 1.0 at the SOL and PFR edges. 
Each returning atom has an energy of 3.5eV:

.. code-block:: ini

   [hermes]
   components = d+, d, sheath_boundary, recycling

   [recycling]
   species = d+   # Comma-separated list of species to recycle

   [d+]
   recycle_as = d         # Species to recycle as

   target_recycle = true  
   target_recycle_multiplier = 0.95 # Recycling fraction
   target_recycle_energy = 3.5   # Energy of recycled particles [eV]

   sol_recycle = true
   sol_recycle_multiplier = 1 # Recycling fraction
   sol_recycle_energy = 3.5   # Energy of recycled particles [eV]

   pfr_recycle = true
   pfr_recycle_multiplier = 1 # Recycling fraction
   pfr_recycle_energy = 3.5   # Energy of recycled particles [eV]

Allowing for fast recycling
^^^^^^^^^^^^^^^

In reality, a fraction of incident ions will undergo specular reflection off the surface and 
preserve a fraction of their energy. In the popular Monte-Carlo neutral code EIRENE, the 
fast recycling fraction and the energy reflection factor are provided by the `TRIM database <https://www.eirene.de/old_eirene/html/surface_data.html>`_
as a function of incident angle, surface material and incident particle energy.
Studies found that sheath acceleration can make the ion angle relatively consistent, e.g. 60 degrees; in (`Jae-Sun Park et al 2021 Nucl. Fusion 61 016021 <https://iopscience.iop.org/article/10.1088/1741-4326/abc1ce>`_).

The recycled heat flux is:

.. math::

   \begin{aligned}
   \Gamma_{E_{n}} &= R \times (R_{f} \alpha_{E} \Gamma_{E_{i}}^{sheath}  + (1 - R_{f}) T_{R} \Gamma_{N_{i}})) \\
   \end{aligned}

Where :math:`R` is the recycle multiplier, :math:`R_{f}` is the fast reflection fraction, :math:`\alpha_{E}` is the energy reflection factor,
:math:`\Gamma_{E_{i}}^{sheath}` is the incident heat flux from the sheath boundary condition, :math:`T_{R}` is the recycle energy and :math:`\Gamma_{N_{i}}` is the incident ion flux.

:math:`R_{f}` and :math:`\alpha_{E}` can be set as in the below example. They can also be set to different values for the SOL and PFR by replacing
the word "target" with either "sol" or "pfr".

.. code-block:: ini

   [d+]
   recycle_as = d         # Species to recycle as

   target_recycle = true  
   target_recycle_multiplier = 0.95 # Recycling fraction
   target_recycle_energy = 3.5   # Energy of recycled particles [eV]
   target_fast_recycle_energy_factor = 0.70
   target_fast_recycle_fraction = 0.80

Neutral pump
^^^^^^^^^^^^^^^

The recycling component also features a neutral pump which is currently implemented for 
the SOL and PFR edges only, and so is not available in 1D. The pump is a region of the wall
which facilitates particle loss by incomplete recycling and neutral absorption. 

The pump requires wall recycling to be enabled on the relevant wall region.

The particle loss rate :math:`\Gamma_{N_{n}}` is the sum of the incident ions that are not recycled and the 
incident neutrals which are not reflected, both of which are controlled by the pump multiplier :math:`M_{p}` 
which is set by the `pump_multiplier` option in the input file. The unrecycled ion flux :math:`\Gamma_{N_{i}}^{unrecycled}` is calculated using the recycling
model and allows for either thermal or fast recycling, but with the difference that the `pump_multiplier` replaces the `recycle_multiplier`. 

.. math::

   \begin{aligned}
   \Gamma_{N_{n}} &= \Gamma_{N_{i}}^{unrecycled} + M_{p} \times \Gamma_{N_{n}}^{incident} \\
   \Gamma_{N_{n}}^{incident} &= N_{n} v_{th} = N_{n} \frac{1}{4} \sqrt{\frac{8 T_{n}}{\pi m_{n}}} \\
   \end{aligned}

Where the thermal velocity formulation is for a static maxwellian in 1D (see Stangeby p.64, eqns 2.21, 2.24) 
and the temperature is in `eV`.

The heat loss rate :math:`\Gamma_{E_{n}}` is calculated as:

.. math::

   \begin{aligned}
   \Gamma_{E_{n}} &= \Gamma_{E_{i}}^{unrecycled}  + M_{p} \times \Gamma_{E_{n}}^{incident} \\
   \Gamma_{E_{n}}^{incident} &= \gamma T_{n} N_{n} v_{th} = 2 T_{n} N_{n} \frac{1}{4} \sqrt{\frac{8 T_{n}}{\pi m_{n}}} \\
   \end{aligned}

Where the incident heat flux is for a static maxwellian in 1D (see Stangeby p.69, eqn 2.30).

The pump will be placed in any cell that
 1. Is the final domain cell before the guard cells
 2. Is on the SOL or PFR edge
 3. Has a `is_pump` value of 1

The field `is_pump` must be created by the user and added to the grid file as a `Field2D`.

Diagnostic variables
^^^^^^^^^^^^^^^
Diagnostic variables for the recycled particle and energy fluxes are provided separately for the targets, the pump as well as the SOL and PFR which are grouped together as `wall`.
as well as the pump. In addition, the field `is_pump` is saved to help in plotting the pump location.


.. doxygenstruct:: Recycling
   :members:
      
.. _binormal_stpm:

Others
~~~~~~~~~~~~~~~

noflow_boundary
^^^^^^^^^^^^^^^

This is a species component which imposes a no-flow boundary condition
on y (parallel) boundaries.

- Zero-gradient boundary conditions are applied to `density`,
  `temperature` and `pressure` fields, if they are set.
- Zero-value boundary conditions are applied to `velocity` and
  `momentum` if they are set.

By default both yup and ydown boundaries are set, but can be turned
off by setting `noflow_lower_y` or `noflow_upper_y` to `false`.

Example: To set no-flow boundary condition on an ion `d+` at the lower
y boundary, with a sheath boundary at the upper y boundary:

.. code-block:: ini

   [hermes]
   components = d+, sheath_boundary

   [d+]
   type = noflow_boundary

   noflow_lower_y = true   # This is the default
   noflow_upper_y = false  # Turn off no-flow at upper y for d+ species

   [sheath_boundary]
   lower_y = false         # Turn off sheath lower boundary for all species
   upper_y = true

Note that currently `noflow_boundary` is set per-species, whereas
`sheath_boundary` is applied to all species. This is because sheath
boundary conditions couple all charged species together, and doesn't
affect neutral species.

The implementation is in `NoFlowBoundary`:

.. doxygenstruct:: NoFlowBoundary
   :members:

.. _neutral_boundary:

neutral_boundary
^^^^^^^^^^^^^^^

Sets Y (sheath/target) boundary conditions on neutral particle
density, temperature and pressure. A no-flow boundary condition
is set on parallel velocity and momentum. It is a species-specific
component and so goes in the list of components for the species
that the boundary condition should be applied to.

Just like ions can undergo fast and thermal recycling, neutrals can undergo fast or thermal 
reflection at the wall. In edge codes using the kinetic neutral code EIRENE, this is typically
controlled by the `TRIM database <https://www.eirene.de/old_eirene/html/surface_data.html>`_.
Hermes-3 features a simpler implementation for a constant, user-set fast reflection fraction :math:`R_{f}`
and energy reflection coefficient :math:`\alpha_{n}` based on the approach in the thesis of D.Power 2023.

The two types of reflection are as follows:

- Fast reflection, where a neutral atom hits the wall and reflects having lost some energy,
- Thermal reflection, where a neutral atom hits the wall, recombines into a molecule, and then
  is assumed to immediately dissociate at the Franck Condon dissociation temperature of 3eV.

They are both implemented as a neutral energy sink calculated
from the cooling heat flux :math:`Q_{cool}`:

.. math::
   \begin{aligned}
   Q_{cool} &= Q_{inc} - Q_{fast_refl} - Q_{th_refl}  \\
   Q_{incident} &= 2n_{n} T_{n} v_{th}^{x}  \\
   Q_{fast} &= 2n_{n} T_{n} v_{th}^{x} (R_{f} \alpha_{n}) \\
   Q_{thermal} &= T_{FC} n_{n} v_{th}^{x} (1 - R_{f}) \\
   v_{th}^{x} &= \frac{1}{4}\sqrt{\frac{8k_{B}T_{n}}{\pi m_{n}}}
   \end{aligned}

Where :math:`Q_{incident}` is the neutral heat flux incident on the wall, :math:`Q_{fast}` is the
returning heat flux from fast reflection, :math:`Q_{thermal}` is the returning heat flux from thermal reflection
and :math:`T_{FC}` is the Franck-Condon dissociation temperature, currently hardcoded to 3eV.
Note that the fast and incident heat flux are both of a Maxwellian distribution, and so their
formula corresponds to the 1 dimensional static Maxwellian heat flux and :math:`v_{th}^{x}` the 
corresponding 1D static Maxwellian thermal velocity (Stangeby p.69).
The thermal heat flux represents a monoenergetic distribution at :math:`T_{n}=T_{FC}` and 
is therefore calculated with a simpler formula.


Since different regions of the tokamak feature different incidence angles and may feature 
different materials, the energy reflection coefficient and the fast reflection fraction 
can be set individually for the target, PFR and SOL walls. The default values are 0.75
for :math:`\alpha_{n}` and 0.8 for :math:`R_{r}` and correspond to approximate values for 
tungsten for incidence angles seen at the target. (Power, 2023)

Here are the options set to their defaults. Note that the SOL and PFR are set to have no
reflection by default so that it is compatible with a model of any dimensionality which has a target.

.. code-block:: ini

   [hermes]
   components = d

   [d]
   type = ... , neutral_boundary

   neutral_boundary_sol = true
   neutral_boundary_pfr = true
   neutral_boundary_upper_y = true
   neutral_boundary_lower_y = true 

   target_energy_refl_factor = 0.75
   sol_energy_refl_factor = 0.75
   pfr_energy_refl_factor = 0.75

   target_fast_refl_fraction = 0.80
   sol_fast_refl_fraction = 0.80
   pfr_fast_refl_fraction = 0.80

.. doxygenstruct:: NeutralBoundary
   :members:

Sources
-------------------

Applying sources using the input file
~~~~~~~~~~~~~~~
The simplest way to implement a source in one of the Hermes-3 equations is through the input file.
This is done by defining an array representing values of the source across the entire domain
using the BOUT++ input file syntax (see `BOUT++ documentation
<https://bout-dev.readthedocs.io/en/latest/user_docs/bout_options.html>`_).

Sources are available for the density, pressure and momentum equations, and are prescribed under 
a header corresponding to the chosen equation and species.

For example, this is how a pressure source is prescribed in the 1D-recycling example. First the domain and grid
are defined using input file functions. This creates a 400 element 1D grid with a length of 30m and an X-point at the 10m mark.
The grid increases in resolution towards the target, with a minimum grid spacing of 0.1 times the average grid spacing:

.. code-block:: ini
   
   [mesh]
   # 1D simulation, use "y" as the dimension along the fieldline
   nx = 1
   ny = 400   # Resolution along field-line
   nz = 1
   length = 30           # Length of the domain in meters
   length_xpt = 10   # Length from midplane to X-point [m] (i.e. this is where the source ends)

   dymin = 0.1  # Minimum grid spacing near target, as fraction of average. Must be > 0 and < 1

   # Parallel grid spacing — grid refinement near the divertor target (which is where the interesting
   # stuff happens)
   dz = 1
   dy = (length / ny) * (1 + (1-dymin)*(1-y/pi))

   # Calculate where the source ends in grid index (i.e. at the X-point)
   source = length_xpt / length
   y_xpt = pi * ( 2 - dymin - sqrt( (2-dymin)^2 - 4*(1-dymin)*source ) ) / (1 - dymin)

And here is how the calculated geometric information is used to prepare a pressure source. The user 
inputs a parallel heatflux in :math:`W/m^2`, or Watts per cross-sectional flux tube area.
This is converted to a pressure flux in :math:`Pa/{m^2s}` by the :math:`2/3` factor, and then
converted to a pressure source in :math:`Pa/{m^3s}` by dividing by the length of the heating region ``mesh:length_xpt``. 
Note that this assumes a constant cross-sectional area, i.e. :math:`dx = dz = J = 1` due to the fact this 
is a 1D simulation. Note that ``dz`` is actually :math:`2 \pi` by default, and must be set to 1 in the input file for
this particular expression to work.
If you are imposing a full B-field profile in your 1D simulation, you will need to account for the fact that :math:`J` is no longer constant.
In order to limit the pressure source to just the region above the X-point, it is multiplied by a Heaviside
function which returns 1 upstream of :math:`y=mesh:y\_xpt` and 0 downstream of it.

.. code-block:: ini

   [Pd+]

   # Initial condition for ion pressure (in terms of hermes:Nnorm * hermes:Tnorm)
   function = 1

   # Input power flux to ions in W/m^2
   powerflux = 2.5e7

   source = (powerflux*2/3 / (mesh:length_xpt))*H(mesh:y_xpt - y)  # Input power as function of y

   [Pe]

   function = `Pd+:function`  # Same as ion pressure initially

   # Input power flux to electrons in W/m^2
   source = `Pd+:source`  # Same as ion pressure source

Applying sources using the grid file
~~~~~~~~~~~~~~~
The input file has limitations, and sometimes it is useful to prepare an arbitrary profile outside of BOUT++
and import it through the grid file. In 2D, this can be done by adding an appropriate Field3D or Field2D to the
grid netCDF file with the sources in the appropriate units.

Time-dependent sources
~~~~~~~~~~~~~~~
Any source can be made time-dependent by adding a flag and providing a prefactor function in the input file.
The already defined source will be multiplied by the prefactor, which is defined by a time-dependent input file function.

Here is the implementation in the 1D-time-dependent-sources example, where the electrons and ions are set to receive 8MW
of mean power flux each with a +/-10% sinusoidal fluctuation of a period of 50us. The density source has a mean of zero and 
oscillates between :math:`-1\times10^{22}` and :math:`1\times10^{22}`, also with a period of 50us.

Note that if you have the density controller enabled, it will work to counteract the imposed density source oscillation.

.. code-block:: ini

   [Nd+]
   function = 5e19 / hermes:Nnorm # Initial conditions
   source_time_dependent = true
   source = 1e22 * H(mesh:y_xpt - y)
   source_prefactor = sin((2/50)*pi*1e6*t)   #  Oscillation between -1 and 1, period 50us

   [Pe]
   function = 0.01
   powerflux = 16e6  # Input power flux in W/m^2
   source = 0.5 * (powerflux*2/3 / (mesh:length_xpt))*H(mesh:y_xpt - y)  # Input power as function of y
   source_time_dependent = true
   source_prefactor = 1 + 0.1 * sin((2/50)*pi*1e6*t)   #  10% fluctuation on on  top of background source, period 50us

   [Pd+]
   function = 0.01
   source = Pe:source
   source_time_dependent = true
   source_prefactor = Pe:source_prefactor