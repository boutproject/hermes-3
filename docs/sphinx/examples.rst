.. _sec-examples:

Examples
========

Hermes-3 contains a large number of example simulations, with only a selection
shown here. Please refer to the ``examples`` directory for more.

As the development of the code has progressed rapidly, it has been difficult
To keep all of the examples up to date and not all are guaranteed to represent
the most optimal setup.

A review is currently in process to select a few of the most useful examples
and keep them up to date regularly, with the remaining examples being updated
in the future.

So far, only 1D example cases received this treatment, resulting in a new,
updated example ``1D-threshold`` and the remaining 1D cases moved to ``examples/1D-tokamak/extra``.
Note that many will likely work with no modifications.

1D Field line
-----------------------

.. _1D-threshold:

1D-threshold
~~~~~~

This is the 
This simulates a similar setup to the `SD1D
<https://github.com/boutproject/SD1D/>`_ code: A 1D domain, with a
source of heat and particles on one side, and a sheath boundary on the
other. Ions recycle into neutrals, which charge exchange and are
ionised.  A difference is that separate ion and electron temperatures
are evolved here.

.. figure:: figs/1d_threshold.*
   :name: 1d_threshold
   :alt:
   :width: 60%

   Evolution of ion and neutral density (blue); ion, electron and
   neutral temperature (red), starting from flat profiles.

Due to the short length-scales near the sheath, the grid is packed
close to the target, by setting the grid spacing to be a linear
function of index:

.. code-block:: ini

   [mesh]
   dy = (length / ny) * (1 + (1-dymin)*(1-y/pi))

where `dymin` is 0.1 here, and sets the smallest grid spacing (at the
target) as a fraction of the average grid spacing.

The components are ion species `d+`, atoms `d`, electrons `e`:

.. code-block:: ini

   [hermes]
   components = (d+, d, e,
              sheath_boundary_simple, collisions, recycling, reactions,
              electron_force_balance, neutral_parallel_diffusion)

The electron velocity is set to the ion by specifying :ref:`zero_current`;
A sheath boundary is included; Collisions are needed to be able to calculate
heat conduction, as well as neutral diffusion rates; Recycling at the targets
provides a source of atoms; :ref:`neutral_parallel_diffusion` simulates cross-field
diffusion in a 1D system. The electron force balance links electron pressure gradient
with the ion momentum equation. Please see the relevant documentation pages about these
components for further information.

The sheath boundary is only imposed on the upper Y boundary:

.. code-block:: ini

   [sheath_boundary]

   lower_y = false
   upper_y = true

The reactions component is a group, which lists the reactions included:

.. code-block:: ini

   [reactions]
   type = (
           d + e -> d+ + 2e,   # Deuterium ionisation
           d+ + e -> d,          # Deuterium recombination
           d + d+ -> d+ + d,   # Charge exchange
          )

To run this example:

.. code-block:: bash

   nice -n 10 ./hermes-3 -d examples/1D-recycling

This should take 5-10 minutes to run. There is a `makeplots.py` script in the
`examples/1D-recycling` directory which will generate plots and a gif animation
(if `ImageMagick <https://imagemagick.org/index.php>`_ is installed).


.. 
   1D-te-ti
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   A fluid is evolved in 1D, imposing quasineutrality and zero net current.
   Both electron and ion pressures are evolved, but there is no exchange
   of energy between them, or heat conduction.

   .. figure:: figs/1d_te_ti.*
      :name: 1d_te_ti
      :alt:
      :width: 60%
      
      Evolution of pressure, starting from a top hat. Input in ``examples/1D-te-ti``.

   To run this example:

   .. code-block:: bash

      ./hermes-3 -d examples/1D-te-ti

   Which takes a few seconds to run on a single core. Then in the
   ``examples/1D-te-ti`` directory run the analysis script

   .. code-block:: bash

      python3 makeplot.py

   That should generate png files and an animated gif if ImageMagick is
   installed (the ``convert`` program). If an error like
   ``ModuleNotFoundError: No module named 'boutdata'`` occurs, then
   install the ``boutdata`` package with ``python3 -m pip install
   boutdata``.

   The model components are ions (i) and electrons (e), and a component
   which uses the force on the electrons to calculate the parallel electric field,
   which transfers the force to the ions.

   .. code-block:: ini

      [hermes]
      components = i, e, electron_force_balance


   The ion density, pressure and momentum equations are evolved:

   .. code-block:: ini

      [i]  # Ions
      type = evolve_density, evolve_pressure, evolve_momentum

   which solves the equations

   .. math::

      \begin{aligned}
      \frac{\partial n_i}{\partial t} =& -\nabla\cdot\left(n_i\mathbf{b}v_{||i}\right) \\
      \frac{\partial p_i}{\partial t} =& -\nabla\cdot\left(p_i\mathbf{b}v_{||i}\right) - \frac{2}{3}p_i\nabla\cdot\left(\mathbf{b}v_{||i}\right) \\
      \frac{\partial}{\partial t}\left(n_iv_{||i}\right) =& -\nabla\cdot\left(n_iv_{||i} \mathbf{b}v_{||i}\right) - \partial_{||}p_i + E
      \end{aligned}

   The electron density is set to the ion density by quasineutrality, the
   parallel velocity is set by a zero current condition, and only the
   electron pressure is evolved.

   .. code-block:: ini

      [e] # Electrons
      type = quasineutral, zero_current, evolve_pressure

   which adds the equations:

   .. math::

      \begin{aligned}
      n_e =& n_i \\
      \frac{\partial p_e}{\partial t} =& -\nabla\cdot\left(p_e\mathbf{b}v_{||e}\right) - \frac{2}{3}p_e\nabla\cdot\left(\mathbf{b}v_{||e}\right)
      \end{aligned}

   The :ref:`zero_current` component sets:

   .. math::

      \begin{aligned}
      E =& -\partial_{||}p_e \\
      v_{||e} =& v_{||i}
      \end{aligned}


2D Drift-plane
--------------

Simulations where the dynamics along the magnetic field is not
included, or only included in a parameterised way as sources or
sinks. The field line direction is then "into the page", and the
domain represents a slice somewhere along the field line, e.g. 
at the midplane.
 These are useful for the study of the basic physics of plasma
"blobs" / filaments, and tokamak edge turbulence. 

.. _Blob2d:

Blob2d
~~~~~~

A seeded plasma filament in 2D. This version is isothermal and cold ion,
so only the electron density and vorticity are evolved. A sheath-connected
closure is used for the parallel current.

.. figure:: figs/blob2d.png
   :name: fig-blob2d
   :alt:
   :scale: 50
   
   Electron density Ne at three times, showing propagation to the right

The model components are

.. code-block:: ini

   [hermes]
   components = e, vorticity, sheath_closure

The electron component consists of two types:

.. code-block:: ini

   [e]  # Electrons
   type = evolve_density, isothermal


The :ref:`evolve_density` component type evolves the electron density `Ne`. This component
has several options, which are set in the same section e.g.

.. code-block:: ini

   poloidal_flows = false  # Y flows due to ExB

and so solves the equation:

.. math::

   \begin{aligned}
   \frac{\partial n_e}{\partial t} =& - \nabla\cdot\left(n_e\mathbf{v}_{E\times B}\right) + \nabla\cdot{\frac{1}{e}\mathbf{j}_{sh}}
   \end{aligned}

The :ref:`isothermal` component type sets the temperature to be a constant, and using
the density then sets the pressure. The constant temperature is also
set in this `[e]` section:

.. code-block:: ini

   temperature = 5  # Temperature in eV

so that the equation solved is

.. math::

   \begin{aligned}
   p_e =& e n_e T_e
   \end{aligned}

where :math:`T_e` is the fixed electron temperature (5eV).

The :ref:`vorticity` component uses the pressure to calculate the diamagnetic current,
so must come after the `e` component. This component then calculates the potential.
Options to control the vorticity component are set in the `[vorticity]` section.

.. math::

   \begin{aligned}
   \frac{\partial \omega}{\partial t} =& - \nabla\cdot\left(\omega\mathbf{v}_{E\times B}\right) + \nabla\left(p_e\nabla\times\frac{\mathbf{b}}{B}\right) + \nabla\cdot\mathbf{j}_{sh} \\
   \nabla\cdot\left(\frac{1}{B^2}\nabla_\perp\phi\right) = \omega
   \end{aligned}

The `sheath_closure` component uses the potential, so must come after :ref:`vorticity`.
Options are also set as

.. code-block:: ini

   [sheath_closure]
   connection_length = 10 # meters

This adds the equation

.. math::

   \begin{aligned}
   \nabla\cdot{\mathbf{j}_{sh}} = \frac{n_e\phi}{L_{||}}
   \end{aligned}

where :math:`L_{||}` is the connection length.

.. _Blob2d-Te-Ti:

Blob2D-Te-Ti
~~~~~~~~~~~~

A seeded plasma filament in 2D. This version evolves both electron and
ion temperatures. A sheath-connected closure is used for the parallel
current.

.. figure:: figs/blob2d-te-ti.png
   :name: fig-blob2d-te-ti
   :alt:
   :scale: 50
   
   Electron density Ne at three times, showing propagation to the right and downwards

The model components are

.. code-block:: ini

   [hermes]
   components = e, h+, vorticity, sheath_closure


The electron component evolves density (saved as `Ne`) and pressure
(`Pe`), and from these the temperature is calculated.

.. code-block:: ini

   [e]
   type = evolve_density, evolve_pressure


The ion component sets the ion density from the electron density, by
using the quasineutrality of the plasma; the ion pressure (`Ph+`) is evolved.

.. code-block:: ini
   
   [h+]
   type = quasineutral, evolve_pressure

The equations this solves are similar to the previous :ref:`Blob2d` case, except
now there are pressure equations for both ions and electrons:

.. math::

   \begin{aligned}
   \frac{\partial n_e}{\partial t} =& - \nabla\cdot\left(n_e\mathbf{v}_{E\times B}\right) + \nabla\cdot{\frac{1}{e}\mathbf{j}_{sh}} \\
   \frac{\partial p_e}{\partial t} =& - \nabla\cdot\left(p_e\mathbf{v}_{E\times B}\right) - \gamma_e p_e c_s \\
   n_{h+} =& n_e \\
   \frac{\partial p_{h+}}{\partial t} =& - \nabla\cdot\left(p_{h+}\mathbf{v}_{E\times B}\right) \\
   \frac{\partial \omega}{\partial t} =& - \nabla\cdot\left(\omega\mathbf{v}_{E\times B}\right) + \nabla\left[\left(p_e + p_{h+}\right)\nabla\times\frac{\mathbf{b}}{B}\right] + \nabla\cdot\mathbf{j}_{sh} \\
   \nabla\cdot\left[\frac{1}{B^2}\nabla_\perp\left(\phi + p_{h+}\right)\right] =& \omega \\
   \nabla\cdot{\mathbf{j}_{sh}} =& \frac{n_e\phi}{L_{||}}
   \end{aligned}

2D-drift-plane-turbulence-te-ti
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A 2D turbulence simulation, similar to the :ref:`Blob2d-Te-Ti` case, but with
extra source and sink terms, so that a statistical steady state of
source-driven turbulence can be reached.

The model components are

.. code-block:: ini

   [hermes]
   components = e, h+, vorticity, sheath_closure


The electron component evolves density (saved as `Ne`) and pressure
(`Pe`), and from these the temperature is calculated.

.. code-block:: ini

   [e]
   type = evolve_density, evolve_pressure


The ion component sets the ion density from the electron density, by
using the quasineutrality of the plasma; the ion pressure (`Ph+`) is evolved.

.. code-block:: ini

   [h+]
   type = quasineutral, evolve_pressure

The sheath closure now specifies that additional sink terms should be added

.. code-block:: ini

    [sheath_closure]
    connection_length = 50 # meters
    potential_offset = 0.0  # Potential at which sheath current is zero
    sinks = true

and radially localised sources are added in the `[Ne]`, `[Pe]`, and `[Ph+]`
sections.

The equations this solves are the same as the previous
:ref:`Blob2d-Te-Ti` case, except wih extra source and sink terms. In
SI units (except temperatures in eV) the equations are:

.. math::

   \begin{aligned}
   p_\mathrm{total} =& \sum_a e n_a T_a \\
   \rho_\mathrm{total} =& \sum_a A_a m_p n_a \\
   c_s =& \sqrt{\frac{p_\mathrm{total}}{\rho_\mathrm{total}}} \\
   \frac{\partial n_e}{\partial t} =& - \nabla\cdot\left(n_e\mathbf{v}_{E\times B}\right) + \nabla\cdot{\frac{1}{e}\mathbf{j}_{sh}} - \frac{n_e c_s}{L_{||}} + S_n \\
   \frac{\partial p_e}{\partial t} =& - \nabla\cdot\left(p_e\mathbf{v}_{E\times B}\right) - \frac{\gamma_e p_e c_s}{L_{||}} + S_{p_e} \\
   n_{h+} =& n_e \\
   \frac{\partial p_{h+}}{\partial t} =& - \nabla\cdot\left(p_{h+}\mathbf{v}_{E\times B}\right) - \frac{\gamma_i p_{h+} c_s}{L_{||}} + S_{p_{h+}} \\
   \frac{\partial \omega}{\partial t} =& - \nabla\cdot\left(\omega\mathbf{v}_{E\times B}\right) + \nabla\cdot\left[\left(p_e + p_{h+}\right)\nabla\times\frac{\mathbf{b}}{B}\right] + \nabla\cdot\mathbf{j}_{sh} \\
   \nabla\cdot\left[\frac{\overline{A}m_p}{B^2}\left(\overline{n}\nabla_\perp\phi + \nabla_\perp p_{h+}\right)\right] =& \omega \\
   \nabla\cdot{\mathbf{j}_{sh}} =& \frac{e n_e \overline{c_s} \phi}{\overline{T} L_{||}} \\
   \mathbf{v}_{E\times B} =& \frac{\mathbf{B}\times\nabla\phi}{B^2}
   \end{aligned}

Where :math:`\overline{T}` and :math:`\overline{n}` are the reference
temperature (units of eV) and density (in units of :math:`m^{-3}`)
used for normalisation. :math:`\overline{c_s} = \sqrt{e\overline{T} /
m_p}` is the reference sound speed, where :math:`m_p` is the proton
mass. The mean ion atomic mass :math:`\overline{A}` is set to 1 here.

These reference values enter into the sheath current
:math:`\mathbf{j}_{sh}` because that is a simplified, linearised form
of the full expression. Likewise the vorticity (:math:`\omega`)
equation used the Boussinesq approximation to simplify the
polarisation current term, leading to constant reference values being
used.

The sheath heat transmission coefficients default to :math:`\gamma_e = 6.5` and
:math:`\gamma_i = 2.0` (:math:`\gamma_i` as suggested in Stangeby's textbook
between equations (2.92) and (2.93)). Note the sinks in may not be correct or
the best choices, especially for cases with multiple ion species; they were
chosen as being simple to implement by John Omotani in May 2022.


2D Axisymmetric SOL
-----------------------

These are transport simulations, where the cross-field transport is given
by diffusion, and fluid-like equations are used for the parallel dynamics
(as in the 1D flux tube cases).

The input settings (in BOUT.inp) are set to read the grid from a file `tokamak.nc`.
This is linked to a default file `compass-36x48.grd.nc`, a COMPASS-like lower single
null tokamak equilibrium. Due to the way that BOUT++ uses communications between
processors to implement branch cuts, these simulations require a multiple of 6 processors.
You don't usually need 6 physical cores to run these cases, if MPI over-subscription
is enabled.

heat-transport
~~~~~~~~~~~~~~

In `examples/tokamak/heat-transport`, this evolves only electron pressure with
a fixed density. It combines cross-field diffusion with parallel heat conduction
and a sheath boundary condition.

To run this simulation with the default inputs requires (at least)
6 processors because it is a single-null tokamak grid.
From the build directory:

.. code-block:: bash

   cd examples/tokamak
   mpirun -np 6 ../../hermes-3 -d heat-transport

That will read the grid from `tokamak.nc`, which by default links to
the `compass-36x48.grd.nc` file.

The components of the model are given in `heat-transport/BOUT.inp`:

.. code-block:: ini

   [hermes]
   components = e, h+, collisions, sheath_boundary_simple

We have two species, electrons and hydrogen ions, and add collisions
between them and a simple sheath boundary condition.

The electrons have the following components to fix the density,
evolve the pressure, and include anomalous cross-field diffusion:

.. code-block:: ini

   [e]
   type = fixed_density, evolve_pressure, anomalous_diffusion

The `fixed_density` takes these options:

.. code-block:: ini

   AA = 1/1836
   charge = -1
   density = 1e18 # Fixed density [m^-3]

so in this simulation the electron density is a uniform and constant value.
If desired, that density can be made a function of space (`x` and `y` coordinates).

The `evolve_pressure` component has thermal conduction enabled, and outputs
extra diagnostics i.e. the temperature `Te`:

.. code-block:: ini

   thermal_conduction = true   # Spitzer parallel heat conduction
   diagnose = true   # Output additional diagnostics

There are other options that can be set to modify the behavior,
such as setting `kappa_limit_alpha` to a value between 0 and 1 to impose
a free-streaming heat flux limit.

Since we're evolving the electron pressure we should set initial and
boundary conditions on `Pe`:

.. code-block:: ini

   [Pe]
   function = 1
   bndry_core = dirichlet(1.0)  # Core boundary high pressure 
   bndry_all = neumann

That sets the pressure initially uniform, to a normalised value of 1,
and fixes the pressure at the core boundary. Other boundaries are set
to zero-gradient (neumann) so there is no cross-field diffusion of heat out of
the outer (SOL or PF) boundaries. Flow of heat through the sheath is
governed by the `sheath_boundary_simple` top-level component.

The hydrogen ions need a density and temperature in order to calculate
the collision frequencies. If the ion temperature is fixed to be the same
as the electron temperature then there is no transfer of energy between
ions and electrons:

.. code-block:: ini

   [h+]
   type = quasineutral, set_temperature

The `quasineutral` component sets the ion density so that there is no net charge
in each cell. In this case that means the hydrogen ion density is set equal to
the electron density. To perform this calculation the component requires that the
ion atomic mass and charge are specified:

.. code-block:: ini

   AA = 1
   charge = 1

The `set_temperature` component sets the ion temperature to the temperature of another
species. The name of that species is given by the `temperature_from` option:

.. code-block:: ini

   temperature_from = e  # Set Th+ = Te

The `collisions` component is described in the manual, and calculates both electron-electron
and electron-ion collisions. These can be disabled if desired, using individual options.
There are also ion-ion, electron-neutral, ion-neutral and neutral-neutral collisions that
are not used here.

The `sheath_boundary_simple` component is a simplified Bohm-Chodura sheath boundary
condition, that allows the sheath heat transmission coefficient to be specified for
electrons and (where relevant) for ions.

The equations solved by this example are:

.. math::

   \begin{aligned}
   \frac{3}{2} \frac{\partial P_e}{\partial t} =& \nabla\cdot\left(\kappa_{e||}\mathbf{b}\mathbf{b}\cdot\nabla T_e\right) + \nabla\cdot\left(n_e\chi\nabla_\perp T_e\right) \\
   \kappa_{e||} =& 3.16 P_e \tau_e / m_e \\
   \tau_e =& 1 / \left(\nu_{ee} + \nu_{ei}\right) \\
   \nu_{ee} =& \frac{2 e^4 n_e \ln\Lambda_{ee}}{3\epsilon_0^2 m_e^2 \left(4\pi e T_e / m_e\right)^{3/2}} \\
   \ln\Lambda_{ee} =& 30.4 - \frac{1}{2}\ln n_e + \frac{5}{4}\ln T_e - \sqrt{10^{-5} + \left(\ln T_e - 2\right)^2 / 16} \\
   \nu_{ei} =& \frac{e^4 n_e \ln\Lambda_{ei}\left(1 + m_e / m_i\right)}{3\epsilon_0^2 m_e^2 \left(2\pi T_e (1/m_e + 1/m_i)\right)^{3/2}} \\
   \ln\Lambda_{ei} =& 31 - \frac{1}{2}\ln n_e + \ln T_e
   \end{aligned}

The calculation of the Coulomb logarithms follows the NRL formulary,
and the above expression is used for temperatures above 10eV. See
the `collisions` manual section for the expressions used in other regimes.

recycling-dthene
~~~~~~~~~~~~~~~~
   
The `recycling-dthene` example includes cross-field diffusion,
parallel flow and heat conduction, collisions between species, sheath
boundary conditions and recycling. It simulates the density, parallel
flow and pressure of the electrons; ion species D+, T+, He+, Ne+; and
neutral species D, T, He, Ne.

.. figure:: figs/pe_nvt_nne_2d.png
   :name: recycling-dthene
   :alt:
   :width: 100%

   Electron pressure, parallel tritium flux, and neon atom density. Simulation
   evolves D, T, He, Ne and electron species, including ions and neutral atoms.

The model components are a list of species, and then collective components
which couple multiple species.

.. code-block:: ini

   [hermes]
   components = (d+, d, t+, t, he+, he, ne+, ne, e,
                 collisions, sheath_boundary, recycling, reactions)

Note that long lists like this can be split across multiple lines by
using parentheses. 
                 
Each ion species has a set of components, to evolve the density,
momentum and pressure. Anomalous diffusion adds diffusion of
particles, momentum and energy. For example deuterium ions contain:

.. code-block:: ini
   
   [d+]
   type = evolve_density, evolve_momentum, evolve_pressure, anomalous_diffusion
   AA = 2
   charge = 1

Atomic reactions are specified as a list:

.. code-block:: ini
   
   [reactions]
   type = (
        d + e -> d+ + 2e,   # Deuterium ionisation
        t + e -> t+ + 2e,   # Tritium ionisation
        he + e -> he+ + 2e, # Helium ionisation
        he+ + e -> he,      # Helium+ recombination
        ne + e -> ne+ + 2e, # Neon ionisation
        ne+ + e -> ne,      # Neon+ recombination
       )
