===============================
Output parameters
===============================

The following configuration parameters describe the simulation process and the
output file format. It is possible to simulate the effluent plume at different
release times and record the output in each case.

|

.. confval:: output.variables

    :type: array of strings
    :default: (all variables)

    Variables to include in the output, all describing properties of the plume
    at various times. Possible alternatives are:

    * ``release_time``: Time of release, relative to simulation start [s]
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

.. confval:: output.resolution

    :type: number
    :units: s

    Time resolution of the ``t`` variable in the output file. This is not the
    same as the internal time step, which is chosen automatically by the
    integration algorithm.

|

.. confval:: output.stagnation

    :type: number
    :units: s

    The maximum value of ``t`` in the output file.
    This parameter determines the length of each simulation period.

|

.. confval:: output.frequency

    :type: number
    :units: s

    Time resolution of the ``release_time`` variable in the output file.
    This parameter determines the time interval between subsequent
    plume simulations.

|

.. confval:: output.stop

    :type: number
    :units: s

    The maximum value of ``release_time`` in the output file.
    This parameter determines when the subsequent plume simulations should be
    stopped.

|

.. confval:: output.csv.file

    :type: string

    Write results to the specified file using the CSV format. Output variables
    are comma-separated. Rows are sorted by ``release_time``, then by ``t``.

|

.. confval:: output.nc.file

    :type: string

    Write results to the specified file using the
    `netCDF4 format <https://unidata.github.io/netcdf4-python/>`_. Output
    variables are structured with ``release_time`` as the first
    dimension and ``t`` as the second dimension.
