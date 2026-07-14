# Temperature controller

This example demonstrates the use of the temperature PI controller. It is almost identical to the 1D-recycling example, except that the power source is controlled via the `temperature_feedback` component.

To use this component, your `[e]` section in `BOUT.inp` should have a section like
```
[e] # Electrons
type = ..., temperature_feedback

...

# What temperature should the controller aim for?
temperature_setpoint = 10.0
# If true, the controller will adjust the power source so the target temperature is at
# temperature_setpoint, and if false (default) it will instead control the upstream
# temperature towards temperature_setpoint
control_target_temperature = true
# If true, force the power source to be >= 0. Usually desirable for the power source.
temperature_source_positive = true
# Check Hermes-3 documentation for how this source works. These are the "proportional"
# and "integral" scalars for the controller.
temperature_controller_p = 1e1
temperature_controller_i = 1e-2

# Rather than having independent controllers for each species, we often want to control
# just one species and set the others in terms of that species. This (somewhat complicated)
# section lets you do just that. In species_for_temperature_feedback, you set which species
# should have the PI-calculated source applied to them (in this case, the electrons and the
# deuterium ions). Then, in scaling_factors_for_temperature_feedback, you can set how much
# of the source you want to give to each species. Here, we're setting Pd = Pe, but you
# could feasibly want to set something like Pd+ = 0.5 * Pe, in which case you'd set 
# >>> scaling_factors_for_temperature_feedback = 1.0, 0.5
# Note: the order of species_for_temperature_feedback must match the order of
#       scaling_factors_for_temperature_feedback!!
species_for_temperature_feedback = e, d+
scaling_factors_for_temperature_feedback = 1.0, 1.0

diagnose = true

...
```

Note that you'll also need to rename `source` to `source_shape` in `[Pe]`.
```
[Pe]

# Input power flux to electrons in W/m^2
function = `Pd+:function`  # Same as ion pressure initially

# Input power flux to electron in W/m^2
powerflux = 2.5e7

source_shape = (powerflux*2/3 / (mesh:length_xpt))*H(mesh:y_xpt - y)  # Input power as function of y
```

To visualise the results, there's a simple `plot_control.py` script which shows how the density, temperatures and sources vary as a function of time. The included `.png` is for the first 100 time-steps.