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
Lee, Joseph H. W., and Chu, Vincent H. (2003):
`Turbulent Jets and Plumes - A Lagrangian Approach <https://doi.org/10.1007/978-1-4615-0407-8>`_.

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

The software can be started from the command line as

.. code-block::

    python -m effluent config.toml

or from within python as

.. code-block:: python

    import effluent
    effluent.run("config.toml")

In both cases, simulation details are specified in the
file ``config.toml``, written in the `TOML file format <https://toml.io/en/>`_.
A samle config file is shown below:

.. literalinclude:: ../../tests/examples/neutral/config.toml
    :language: toml

Other input and output formats are also available. See the rest of the
documentation for details.


Documentation
=============
.. toctree::
    :maxdepth: 2

    algorithm
    config
    examples
