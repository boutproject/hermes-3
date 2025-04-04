.. _sec-feedback-control:

Feedback control
===============================

.. _upstream_density_feedback:

Upstream density controller
~~~~~~~~~~~~~~~~~~~~~~~~~

This is intended for 1D simulations, where the density at :math:`y=0` is set
by adjusting an input source. This component uses a PI controller method
to scale the density source up and down, to maintain the specified upstream
density. 
The source, e.g. ``Sd+_feedback``, is calculated as a product of the control signal ``density_source_multiplier``, 
and the array ``density_source_shape`` which defines the source region.
The signal is non-dimensional and the controller depends on the value of ``density_source_shape`` to have a good initial guess of the source.
It should be set to a reasonable value in the units of ``[m-3s-1]``. 
A good reasonable value is the expected steady state domain particle loss (for example due to unrecycled ions at the target).


For example:

.. code-block:: ini

   [d+]
   type = ..., upstream_density_feedback

   density_upstream = 1e19      # Density in m^-3
   density_controller_p = 1e-2  # Feedback controller proportional (p) parameter
   density_controller_i = 1e-3  # Feedback controller integral (i) parameter

   [Nd+]
   source_shape = h(pi - y) * 1e20  # Source shape

There are two additional settings which can make the controller more robust without excessive tuning:

``density_source_positive`` ensures the controller never takes particles away, which can prevent oscillatory
behaviour. Note that this requires some other domain particle sink to ensure control, or else the particle count can never reduce.

``density_integral_positive`` This makes sure the integral component only adds particles. 
The integral component takes a long time to change value, which can result in large overshoots if the initial guess was too small.
This setting mitigates this by disabling the integral term if the density is above the desired value.

Notes:
   - The example cases have their PI parameters tuned properly without the need of the above two settings.
   - Under certain conditions, the use of the PI controller can make the upstream density enter a very small oscillation (~0.05% of upstream value).
   - There is a separate `source` setting that includes a fixed (non varying) density source.

The implementation is in the `UpstreamDensityFeedback` class:

.. doxygenstruct:: UpstreamDensityFeedback
   :members:

Detachment front position controller
~~~~~~~~~~~~~~~~~~~~~~~~~

*N.b. currently this code has been developed for and tested only in 1D
and with a single processor.*

Overview
--------

The ``DetachmentController`` is designed to adjust plasma parameters
dynamically to maintain the detachment front at a desired position. It
uses a proportional-integral-derivative (PID) controller approach,
customized with specific logic for detachment control in plasma
simulations. The controller modifies the source terms for energy or
particle density based on the difference between the current and desired
detachment front locations.

Key components:

-  **Constructor**: Initializes the controller with simulation options,
   including setpoints, control parameters (gain, integral time,
   derivative time), actuator type (power or particles), and other
   relevant settings.

-  **transform(Options& state)**: Main method called at each simulation
   timestep to calculate the control action. It performs several steps:

   1. Computes the current location of the detachment front based on
      plasma densities.
   2. Calculates the error between the desired and current detachment
      front locations.
   3. Updates the control signal based on the PID algorithm.
   4. Applies the control signal to adjust the source terms accordingly.

-  **outputVars(Options& state)**: Outputs diagnostic variables related
   to the control process, such as the detachment front location,
   control signal, and PID terms.

-  **restartVars(Options& state)**: Manages state variables for
   restarting the simulation from saved states, ensuring the controller
   can smoothly continue from a previous point, without perturbations
   due to reinitialising internal variables.

Usage
-----

To use the ``DetachmentController`` in a BOUT++ simulation, you must
provide a setpoint for the detachment front location, choose the
actuator type, and configure the PID controller parameters (gain,
integral time, derivative time). The input for a minimum working example
would be something like

::

   [hermes]
   components = (d+, d, e, detachment_controller, ...)

   [detachment_controller]
   ; Desired location of the detachment front, in metres from the target
   detachment_front_setpoint = 1.0
   ; What is the main neutral species?
   neutral_species = "d+"
   ; What should we adjust to control the detachment front?
   actuator = "power"
   ; What species should we look at for `source_shape`
   species_for_source_shape = e
   ; What species should we control, and how much control should apply
   ; to each species?
   species_list = e, d+
   scaling_factors_list = 1.0, 1.0

   ; What value to use for Kp? See the rest of the documentation
   controller_gain = 1.0

   [Pe]
   function = 1.0
   source_shape = (power_flux*2/3 / (mesh:length_xpt))*H(mesh:y_xpt - y)

   [Pd+]
   function = Pe:function
   ; No Pd+ source or source_shape

Top-level overview of the algorithm
-----------------------------------

-  **Detachment Front Location Calculation**: Determines the point where
   the neutral density exceeds the electron density, which is taken as
   the position of the detachment front.

-  **PID Controller Implementation**: Adjusts the control signal based
   on the error between the current and target detachment front
   locations. Supports both velocity and position form PID control.

-  **Actuator Type**: Defines whether the control action adjusts power
   or particle sources in the simulation.

-  **Control Signal Application**: Modifies the source terms for
   specified plasma species based on the control signal, influencing the
   plasma dynamics to achieve the desired detachment front location.

Fine details
-----------

Finding the detachment front
^^^^^^^^^^^^^^^

We define the detachment front as the point where the neutral density
(``N_neutral_species``) becomes larger than the electron density
(``Ne``). There’s lots of other (better) definitions, but this one is
the easiest to calculate. For each timestep, the controller starts at
the divertor target and checks whether a cell has :math:`n_N > n_e`. The
first time it finds a cell fulfilling that condition, it calculates to
first order the position where the cross-over would occur, by comparing
the last point where :math:`n_N > n_e` (point 1) with the first point
where :math:`n_e > n_N` (point 2);

.. math::


   x_D = \frac{(x_1 \times n_{e,2} - n_{e,1} \times x_2) - (x_1 \times n_{n,2} - n_{n,1}*x_2)}{(n_{n,1} - n_{n,2}) - (n_{e,1} - n_{e,2})}

Note that, since the :math:`y` values aren’t stored directly, we compute
it from :math:`dy` while iterating along cells.

Calculating the feedback response via PID control
^^^^^^^^^^^^^^^

Once we have the detachment position, we need to control its position.
Firstly, we define the error :math:`e = x_S - x_D` where :math:`x_S` is
our setpoint and :math:`x_D` is the actual position of the detachment
front. We also fetch the current time :math:`t`.

If starting a new simulation, we don’t actually want to start
controlling the detachment front straight away, but instead wait for
``settling_time`` (in seconds). Until then, we just set
``control = initial_control``.

Once we’re through the ``settling_time``, we only apply control if the
time has changed by some ``min_time_for_change`` and if the absolute
change of the error is at least ``min_error_for_change``. Usually we set
``min_time_for_change`` to be a small, but non-zero number, which
prevents us from dividing by zero when evaluating the derivative of the
error. Often the ``min_error_for_change`` is zero, although you can
experiment with this and see if it gives you improved performance of the
controller. Once these conditions are met, we update the control value.

There are two options for the control algorithm to use — the ‘position’
and ‘velocity’ forms of the PID. These are, respectively,

.. math::


   C =  C_0 + s\times K_c\times \left[e + \frac{1}{\tau_I} \int e\cdot  dt +\tau_D \frac{d e}{d t} \right]

and

.. math::


   C = C_{prev} + s\times K_c\times \left[\Delta e + \frac{\Delta t}{\tau_I} e +\tau_D \Delta \left(\frac{d e}{d t}\right) \right]

where \* :math:`C` is the control value applied to the source \*
:math:`s` is the response sign (:math:`+1` for particles or :math:`-1`
for power) \* :math:`K_c` is the controller gain \* :math:`\tau_I` is
the integral time \* :math:`\tau_D` is the derivative time \*
:math:`\int e\cdot dt` is the error integral \* :math:`\frac{d e}{d t}`
is the error derivative \* :math:`C_0` is the control offset (usually
zero, unless you’re running in proportional-only mode and need a finite
control when the error is zero) \* :math:`C_{prev}` is the previous
value of :math:`C`

In the position form, we directly update the control value, while in the
velocity form we calculate the change of the control value and add it to
the previous control value. These two approaches should be largely
equivalent, although the velocity form has the advantage of avoiding
integral windup when the control is at either its min or max value.
Conversely, with the position form we can apply other anti-windup
schemes which explicitly change the error integral. A particularly
useful anti-windup method (only available in position form) is
``reset_integral_on_first_crossing``, where the error integral is reset
when the detachment front first reaches the setpoint.

In either scheme, the control value is bounded by
``minval_for_source_multiplier`` and ``maxval_for_source_multiplier``,
which can be used either to ensure that the system stays within
physical, engineering or numerical-stability constraints.

Understanding PID control and tuning the coefficients
^^^^^^^^^^^^^^^

See
`en.wikipedia.org/wiki/Proportional-integral-derivative_controller <https://en.wikipedia.org/wiki/Proportional%E2%80%93integral%E2%80%93derivative_controller>`__.
Most relevant sections are

-  `Controller
   theory <https://en.wikipedia.org/wiki/Proportional%E2%80%93integral%E2%80%93derivative_controller#Controller_theory>`__

   -  `Proportional
      term <https://en.wikipedia.org/wiki/Proportional%E2%80%93integral%E2%80%93derivative_controller#Proportional_term>`__
   -  `Integral
      term <https://en.wikipedia.org/wiki/Proportional%E2%80%93integral%E2%80%93derivative_controller#Integral_term>`__
   -  `Derivative
      term <https://en.wikipedia.org/wiki/Proportional%E2%80%93integral%E2%80%93derivative_controller#Derivative_term>`__

-  `Loop
   tuning <https://en.wikipedia.org/wiki/Proportional%E2%80%93integral%E2%80%93derivative_controller#Loop_tuning>`__
-  `Ziegler-Nichols
   method <https://en.wikipedia.org/wiki/Ziegler%E2%80%93Nichols_method>`__

Input Parameters
^^^^^^^^^^^^^^^

-  **``detachment_front_setpoint``**: The desired position of the
   detachment front from the divertor target, measured in meters from
   the divertor target. It represents the target location where the
   control system aims to maintain the front.

-  **``velocity_form``**: A boolean flag indicating whether to use the
   velocity form (if ``true``) or the position form (if ``false``) of
   the PID controller.

-  **``min_time_for_change``**: The minimum time interval, in seconds,
   before the control signal can be updated. This parameter prevents
   divide-by-zero errors when evaluating the error derivative.

-  **``min_error_for_change``**: The minimum change in error required
   before the control signal is updated. This term is mostly left for
   experimenting.

-  **``minval_for_source_multiplier``**: The minimum value that the
   control signal (source multiplier) can take.

-  **``maxval_for_source_multiplier``**: The maximum value for the
   control signal (source multiplier).

-  **``species_for_source_shape``**: Specifies the plasma species from
   which to select the source shape. The feedback source is the product
   of the control signal (source multiplier) and the source shape. If
   ``actuator='power'`` this will set ``source_shape=Ps::source_shape``
   where ``s=species_for_source_shape``, and if ``actuator='particles'``
   this will set ``source_shape=Ns::source_shape``.

-  **``neutral_species``**: Indicates the main neutral species in the
   plasma. This species is used to determine the location of the
   detachment front.

-  **``actuator``**: Defines the actuator to be adjusted to control the
   detachment front position. Options include ‘power’ for energy sources
   or ‘particles’ for particle density sources.

-  **``initial_control``**: The initial value for the source multiplier
   (control signal) at the start of the simulation.

-  **``control_offset``**: The expected control value when the error
   equals zero. It serves as a baseline around which the control signal
   is modulated. Only used in position form.

-  **``settling_time``**: The time allowed for the system to settle
   before activating certain control terms, measured in seconds. It
   delays the start of control actions to ensure initial transients have
   subsided.

-  **``ignore_restart``**: A flag to ignore the restart file, mainly
   useful for development purposes. It forces the controller to
   initialize from the provided settings rather than loading a previous
   state.

-  **``reset_integral_on_first_crossing``**: Resets the error integral
   to zero when the detachment front first reaches the desired position.
   This feature can help prevent integral wind-up and improve control
   stability.

-  **``controller_gain``** (Kc): The proportional gain of the PID
   controller. It determines the strength of the response to the error.

-  **``integral_time``**: The integral time of the PID controller, which
   influences the rate at which the integral term accumulates error over
   time.

-  **``derivative_time``**: The derivative time of the PID controller,
   affecting how strongly the controller reacts to the rate of change of
   the error.

-  **``buffer_size``**: The number of points to store for calculating
   derivatives. It determines the size of the window over which the
   derivative of the error is computed. Usually a value of around 3 to 5
   seems to filter out the worst of the noise due to small time-steps.
   Very large values will potentially introduce a destabilising lag in
   the derivative response (and also will eventually start to impact
   computational performance).

-  **``species_list``**: A comma-separated list of species to which the
   PI-controlled source will be applied.

-  **``scaling_factors_list``**: A comma-separated list of scaling
   factors corresponding to each species listed in ``species_list``.
   These factors adjust the magnitude of the control action applied to
   each species.

-  **``diagnose``**: Enables the output of additional diagnostic
   variables related to the control process if set to ``true``.

-  **``debug``**: Controls the level of debugging information printed to
   the screen. A value of ``0`` disables debugging output, ``1`` enables
   basic output, and ``2`` provides extensive debugging information.
