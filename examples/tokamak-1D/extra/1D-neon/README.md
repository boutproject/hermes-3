1D scrape-off layer with neon
=============================

This is a more complex version of the `1D-recycling` example. It adds
neon species with multiple charge states, thermal force, ionisation
and recombination between states.

*Restart from 1D-hydrogen* : Use the restart file from the 1D-hydrogen
example as the initial condition for this simulation.

*Solver settings* : The solver settings require PETSc compiled with
STRUMPACK . Solving the system of equations with multiple species
is difficult and requires some tuning. See tips and recipes
here: https://github.com/mikekryjak/hermes-perftest/
