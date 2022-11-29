================================
Neutrally buoyant plume
================================

The centerline of a neutrally buoyant plume stays at the same depth, but the
concentration and momentum is gradually declining as the plume mixes with the
surrounding waters.

To run the example, create input files as described below, and run
:code:`python -m effluent config.toml` from the command line. Alternatively,
execute the statement :code:`effluent.run("config.toml")` from within python.

Input files: **config.toml**, **pipe.csv**, **ambient.csv**

Output file: **out.csv**

Contents of **config.toml**:

.. literalinclude:: ../../../tests/examples/neutral/config.toml
    :language: toml

Contents of **pipe.csv**:

.. literalinclude:: ../../../tests/examples/neutral/pipe.csv

Contents of **ambient.csv**:

.. literalinclude:: ../../../tests/examples/neutral/ambient.csv

Contents of **out.csv**:

.. literalinclude:: ../../../tests/examples/neutral/expected.csv
