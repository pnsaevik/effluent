.. effluent documentation master file, created by
   sphinx-quickstart on Thu Feb 25 10:45:33 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==========================================
Effluent: Pipe discharge model
==========================================


What is effluent?
===================

``effluent`` is a python package for simulating the dispersion of effluent
discharges from wastewater pipes. The underlying model is based on
Lee, Joseph H. W., and Chu, Vincent H. (2003):
`Turbulent Jets and Plumes - A Lagrangian Approach
<https://doi.org/10.1007/978-1-4615-0407-8>`_.

The package is mainly intended for research purposes, and does not contain
plotting or statistics functionality. It is expected that
these analyses are conducted in post-processing stages using other packages.

.. _citation:

Citation
========

If you use the software in a publication or report, please cite it as follows:

Sævik, P. N., (2023). *Effluent: A Python package for modelling effluent discharge*.
Journal of Open Source Software, **8**\ (89),
`<https://doi.org/10.21105/joss.05554>`_

Consider also citing the work which this package is based on:

Lee, J. H. W., and Chu, V. H. (2003).
*Turbulent Jets and Plumes - A Lagrangian Approach*. Boston, MA: Springer US.
`<https://doi.org/10.1007/978-1-4615-0407-8>`_.


Installation
============

The package is installed using pip:

::

  pip install effluent

Usage
=====

The software can be started from the command line as

.. code-block::

    effluent config.toml

or from within python as

.. code-block:: python

    import effluent
    effluent.run("config.toml")

In both cases, simulation details are specified in the
file ``config.toml``, written in the `TOML file format <https://toml.io/en/>`_.
The :ref:`examples_page` section includes many examples of valid config files,
for various types of problems.


Documentation
=============
.. toctree::
    :maxdepth: 2

    algorithm
    config
    examples
    autoapi/index
    contributing
    credits
