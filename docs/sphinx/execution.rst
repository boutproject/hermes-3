.. _sec-execution:

Executing simulations
===========


Running Hermes-3 is done in the same way as any BOUT++ simulation.
Please refer to the relevant `BOUT++ documentation 
<https://bout-dev.readthedocs.io/en/stable/user_docs/running_bout.html>`_
for details, or see below for a simple example.

This command runs the simulation "1D-threshold" using 4 processors:

.. code-block:: ini

  mpirun -np 4 hermes-3/build/hermes-3 -d 1D-threshold

Simulations can be restarted from any previous result:

.. code-block:: ini

  mpirun -np 4 hermes-3/build/hermes-3 -d 1D-threshold restart

And restarted while preserving previous results:

.. code-block:: ini

  mpirun -np 4 hermes-3/build/hermes-3 -d 1D-threshold restart append

Choosing the number of processors
~~~~~~~~~~~

When generating a grid analytically, e.g. in 1D or 2D slab models, you can choose 
an arbitrary number of processors. When using a tokamak grid, the parallel domain
decomposition is hardcoded, and it is not possible to run in serial. Since each 
rank must contain the same number of grid cells, the simulation can only run
on a certain number of processors. The way to calculate this number is rather
obscure but covered in the `BOUT++ documentation 
<https://bout-dev.readthedocs.io/en/stable/user_docs/input_grids.html#advanced>`_. 

Output files
-----------

BOUT++ generates four types of output files:

* ``BOUT.settings``: this will capture all of the options used in the simulation,
  including ones set in BOUT.inp as well as defaults. It is generated only 
  when the simulation completes, and will not be generated if interrupted.

* ``BOUT.log.*``: This captures the entire console output of the simulation, 
  capturing all inputs, any warnings or errors, the version of BOUT++ and
  Hermes-3 used as well as the iteration output. Each process generates
  its own log file, which is named according to the process number.

* ``BOUT.dmp.*``: The dump files are netCDF files containing all of the variables saved
  by the simulation. There is one per process, just like how the log files work.

* ``BOUT.restart.*``: These files are required to restart the simulation, 
  and are generally smaller than the dump files.

.. _sec-execution-squashing:

Squashing output
~~~~~~~~~~~

The multiple dump files can be "squashed" into just one, reducing file size and IO
overhead by using `boutdata.squashoutput
<https://github.com/boutproject/boutdata/blob/
0aaef417af092882ac295c4d84e4532e4a10e01f/src/boutdata/squashoutput.py#L16>`_.

Manipulating restart files
~~~~~~~~~~~

Restart files can be created from an arbitrary time slice of the simulation using 
`boutdata.restart.create <https://github.com/boutproject/
boutdata/blob/0aaef417af092882ac295c4d84e4532e4a10e01f/src/boutdata/restart.py#L459>`_.

