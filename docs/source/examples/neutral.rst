================================
Neutrally buoyant plume
================================

The centerline of a neutrally buoyant plume stays at the same depth, but the
concentration and momentum is gradually declining as the plume mixes with the
surrounding waters.

To run the example, create input files as described below, and run
:code:`python -m effluent config.yaml` from the command line. Alternatively,
execute the statement :code:`effluent.run("config.yaml")` from within python.

Input files: **config.yaml**, **pipe.csv**, **ambient.csv**

Output file: **out.csv**

Contents of **config.yaml**:

.. literalinclude:: ../../../tests/examples/neutral/config.yaml
    :language: yaml

Contents of **pipe.csv**:

.. literalinclude:: ../../../tests/examples/neutral/pipe.csv

Contents of **ambient.csv**:

.. literalinclude:: ../../../tests/examples/neutral/ambient.csv

Contents of **out.csv**:

.. literalinclude:: ../../../tests/examples/neutral/expected.csv
