===============================
Pipe parameters
===============================

Pipe parameters determine the physical properties of
the outlet pipe and the effluent water masses. The
properties may vary with time.

It is required to supply either a complete set of explicit parameters
(:confval:`pipe.time`, :confval:`pipe.flow`, :confval:`pipe.dens`,
:confval:`pipe.diam`, :confval:`pipe.depth` and :confval:`pipe.decline`) or the
name of an external file containing the parameters (:confval:`pipe.csv.file`
or :confval:`pipe.nc.file`).

|

.. confval:: pipe.time

    :type: array
    :units: s

    Time since simulation start.

|

.. confval:: pipe.flow

    :type: array
    :units: m³ / s

    Pipe flow rate.
    One entry for each :confval:`pipe.time` value.

|

.. confval:: pipe.dens

    :type: array
    :units: kg / m³

    Mass density of effluent wastewater.
    One entry for each :confval:`pipe.time` value.

|

.. confval:: pipe.diam

    :type: array
    :units: m

    Diameter of pipe outlet.
    One entry for each :confval:`pipe.time` value.

|

.. confval:: pipe.depth

    :type: array
    :units: m

    Depth of pipe outlet.
    One entry for each :confval:`pipe.time` value.

|

.. confval:: pipe.decline

    :type: array
    :units: degrees

    Direction of pipe outlet (positive is downwards).
    One entry for each :confval:`pipe.time` value.

|

.. confval:: pipe.csv.file

   :type: string

   Read pipe parameters from the specified text file. The file must have one
   column (with header) for each pipe parameter. Columns must be
   comma-separated. Lines starting with ``#`` are treated as comments, and
   whitespace is ignored.

|

.. confval:: pipe.nc.file

   :type: string

   Read pipe parameters from the specified
   `netCDF4 file <https://unidata.github.io/netcdf4-python/>`_.
   The file must have one variable for each pipe parameter, indexed by the time
   coordinate.
