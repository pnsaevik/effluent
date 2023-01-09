=======================================
Pure horizontal jet
=======================================

A pure horizontal jet is the simplest test case, which is also investigated in
numerous experimental studies. The ``config.toml`` file looks like this:

.. literalinclude:: ../../../tests/examples/jet/config.toml
    :language: toml

To run the example, run :code:`python -m effluent config.toml` from the command
line. Alternatively, execute the statement :code:`effluent.run("config.toml")`
from within python.


.. plot::
    :context:
    :nofigs:

    import shutil

    src = "../../../tests/examples/jet/expected.csv"
    dst = "out.csv"
    shutil.copyfile(src, dst)

    src = "../../../tests/examples/jet/config.toml"
    dst = "config.toml"
    shutil.copyfile(src, dst)

.. plot::
    :context:

    import matplotlib.pyplot as plt
    import pandas as pd

    df = pd.read_csv("out.csv")

    diam = 0.5
    flow = 0.2
    area = 3.14 * diam * diam * 0.25
    u0 = flow / area

    dist = 6.2 + df.x.values / 0.5
    u = df.u.values / u0
    plt.loglog(dist, u)



.. plot::
    :context:
    :nofigs:

    from pathlib import Path

    Path("config.toml").unlink()
    Path("out.csv").unlink()
