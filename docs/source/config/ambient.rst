===============================
Ambient parameters
===============================

Ambient parameters determine the physical properties of
the ambient water masses. The properties may with time
and depth.

It is required to supply either a complete set of explicit parameters
(:confval:`ambient.time`, :confval:`ambient.depth`, :confval:`ambient.coflow`,
:confval:`ambient.crossflow` and :confval:`ambient.dens`) or the name of an
external file containing the parameters (:confval:`ambient.csv.file`,
:confval:`ambient.nc.file` or :confval:`ambient.roms.file`).

|

.. confval:: ambient.time

    :type: array
    :units: date

    Time and date associated with the pipe data.

|

.. confval:: ambient.depth

    :type: array
    :units: m

    Depth below sea surface.

|

.. confval:: ambient.coflow

    :type: array of arrays
    :units: m / s

    Horizontal velocity in co-flow direction, i.e., in the direction parallel
    to the pipe.
    There should be one array for each :confval:`ambient.time` value. Each
    array should have the same number of elements as :confval:`ambient.depth`.

|

.. confval:: ambient.crossflow

    :type: array of arrays
    :units: m / s

    Horizontal velocity in crossflow direction. When facing in the same
    direction as the pipe outlet, the positive crossflow direction is to the
    right.
    There should be one array for each :confval:`ambient.time` value. Each
    array should have the same number of elements as :confval:`ambient.depth`.

|

.. confval:: ambient.dens

    :type: array of arrays
    :units: kg / m³

    Mass density of ambient water masses.
    There should be one array for each :confval:`ambient.time` value. Each
    array should have the same number of elements as :confval:`ambient.depth`.

|

.. confval:: ambient.csv.file

    :type: string

    Read ambient parameters from the specified text file. The file must have
    one column (with header) for each ambient parameter. Columns must be
    comma-separated. Lines starting with ``#`` are treated as comments, and
    whitespace is ignored.

|

.. confval:: ambient.nc.file

    :type: string

    Read ambient parameters from the specified
    `netCDF4 file <https://unidata.github.io/netcdf4-python/>`_.
    The file must have one variable for each pipe parameter. Each of the
    two-dimensional variables should have time as its first dimension and depth
    as its second dimension.

|

.. confval:: ambient.roms.file

    :type: string

    Read ambient parameters from the ocean model
    `ROMS <https://www.myroms.org/>`_. If wildcards are given, the matching
    files are assumed to be ordered sequentially in time. The software reads
    the fields ``salt`` and ``temp``, and computes seawater density from them.

    Reading ROMS files requires
    `dask <https://docs.xarray.dev/en/stable/dask.html>`_ to be installed.

    Time must be indexed by the ``ocean_time`` coordinate. Horizontal
    coordinates ``lat_rho`` and ``lon_rho`` must be present in the first file.
    Depth coordinates ``h``, ``Cs_r`` and vertical parameters ``hc``,
    ``Vtransform`` must also be present.

|

.. confval:: ambient.roms.latitude

    :type: number

    Latitude of the relevant data points from :confval:`ambient.roms.file`.

|

.. confval:: ambient.roms.longitude

    :type: number

    Longitude of the relevant data points from :confval:`ambient.roms.file`.

|

.. confval:: ambient.roms.azimuth

    :type: number
    :units: degrees

    Azimuthal direction of the co-flow direction (i.e., the direction of the
    pipe outlet). North is 0 and east is 90.
