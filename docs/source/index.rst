.. effluent documentation master file, created by
   sphinx-quickstart on Thu Feb 25 10:45:33 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==========================================
Effluent: Pipe discharge model
==========================================


What is effluent?
===================

Effluent is a python package for simulating the dispersion of effluent
discharges from wastewater pipes. The underlying model is based on
Lee, Joseph H. W., and Chu, Vincent H.: *Turbulent Jets and Plumes*.
Boston, MA: Springer US, 2003.
`doi:10.1007/978-1-4615-0407-8 <https://doi.org/10.1007/978-1-4615-0407-8>`_.

The package is mainly intended for research purposes, and does not contain
any convenience plotting or statistics functionality. It is expected that
these analyses are conducted in post-processing stages using other packages.


Installation
============

The package is installed using pip:

::

  pip install effluent

Usage
=====

The software is invoked from the command line as

::

  python -m effluent config.yaml

where ``config.yaml`` is the configuration file specifying the simulation
setup. A samle config file is shown below:

.. code-block:: yaml

  # Characteristics of the pipe and the effluent flow
  pipe:
    times: [0, 1200]  # Time since simulation start [s]

    flow: [1, 1]    # [m^3/s]
    temp: [10, 10]  # [degrees celcius]
    salt: [10, 10]  # [psu]

    diam: 0.5   # [m]
    depth: 50   # [m]
    decline: 0  # Direction of outlet (positive is downwards) [degrees]

  # Characteristics of the ambient water masses
  ambient:
    times: [0, 3600]      # Time since simulation start [s]
    depths: [0, 10, 50]   # Depth levels [m]

    # Velocity in co-flow, cross-flow (to the left) and in the upwards
    # directions [m/s]
    coflow: [[.1, .1, .1], [.1, .1, .1]]
    crossflow: [[.2, .2, .2], [.2, .2, .2]]
    upflow: [[0, 0, 0], [0, 0, 0]]

    temp: [[10, 10, 10], [10, 10, 10]]  # [degrees celcius]
    salt: [[34, 34, 34], [34, 34, 34]]  # [kg/m^3]

  # Output file format
  output:
    file: out.nc
    frequency: 3600  # Output frequency [s]

Other input and output formats are also available. See the rest of the
documentation for details.

Documentation
=============
.. toctree::
   :maxdepth: 2

   algo
   input
   output
