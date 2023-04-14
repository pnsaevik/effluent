===============================
Output parameters
===============================

The output file contains a set of trajectories of the effluent jet, one for
each release time. The first release time is
:confval:`output.release.start`, the last release time is
:confval:`output.release.stop`, and the number of intermediate trajectories is
given by :confval:`output.release.step`.

Each trajectory is a collection of points, one for every output time step.
The first trajectory time step is :confval:`output.trajectory.start`, the
last one is :confval:`output.trajectory.stop`, and the number of intermediate
trajectory points is determined by :confval:`output.trajectory.step`. The
number of actual internal time steps is determined by
:doc:`solver parameters </config/solver>`.

Output is either written to a text file (:confval:`output.csv.file`) or a
netCDF file (:confval:`output.nc.file`). The option :confval:`output.variables`
can be used to limit the number of variables written.

|

.. confval:: output.variables

    :type: array of strings
    :default: (all variables)

    Variables to include in the output. Possible alternatives are:

    * ``release_time``: Time of release [date]
    * ``t``: Time since release from pipe outlet [s]
    * ``x``: Centerline *x* position (co-flow direction) [m]
    * ``y``: Centerline *y* position (crossflow direction) [m]
    * ``z``: Centerline *z* position (depth below surface) [m]
    * ``u``: Mean velocity in *x* direction [s]
    * ``v``: Mean velocity in *y* direction [s]
    * ``w``: Mean velocity in *z* direction [s]
    * ``density``: Mean mass density [kg/m³]
    * ``radius``: Radius of top-hat profile [m]

|

.. confval:: output.release.start

    :type: date
    :default: 1970-01-01

    Date and time of the first simulated trajectory.

|

.. confval:: output.release.stop

    :type: date
    :default: 1970-01-01

    Date and time of the last simulated trajectory.

|

.. confval:: output.release.step

    :type: number
    :units: s
    :default: 86400

    Time between each simulated trajectory.

|

.. confval:: output.trajectory.start

    :type: number
    :units: s
    :default: 0

    The first trajectory point (first value of ``t``) written to the output
    file.

|

.. confval:: output.trajectory.stop

    :type: number
    :units: s
    :default: 60

    The last trajectory point (last value of ``t``) written to the output
    file.

|

.. confval:: output.trajectory.step

    :type: number
    :units: s
    :default: 1

    The time between trajectory points (i.e., time between ``t`` values)
    written to the output file. This is not the
    same as the internal time step, which is chosen automatically by the
    integration algorithm.

|

.. confval:: output.csv.file

    :type: string

    Write results to the specified comma-delimited text file.
    Rows are sorted by ``release_time``, then by ``t``.

|

.. confval:: output.csv.float_format

    :type: string
    :default: "%.10g"

    Format and precision of floats written to file. Passed directly to
    `pandas.DataFrame.to_csv <https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_csv.html>`_.

|

.. confval:: output.nc.file

    :type: string

    Write results to the specified file using the
    `netCDF4 format <https://unidata.github.io/netcdf4-python/>`_. Output
    variables are structured with ``release_time`` as the first
    dimension and ``t`` as the second dimension.
