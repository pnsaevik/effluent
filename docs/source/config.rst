==============
Configuration
==============

The software can be started from the command line as

.. code-block::

    python -m effluent config.toml

or from within python as

.. code-block:: python

    import effluent
    effluent.run("config.toml")

In both cases, simulation details are specified in the
file ``config.toml``, written in the `TOML file format <https://toml.io/en/>`_.
Here we describe the different options available:

.. toctree::
    :maxdepth: 2

    config/pipe
    config/ambient
    config/output
    config/model
