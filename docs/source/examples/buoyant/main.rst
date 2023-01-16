=======================================
Buoyant jet in stagnant ambient
=======================================

.. plot::
    :context:

    plt.clf()
    plt.close('all')

|

In this example, a horizontal freshwater jet enters an environment of
stagnant, salty water. The ``config.toml`` file looks like this:

.. literalinclude:: config.toml
    :language: toml

|

The contents of the output file ``out.csv`` is

.. literalinclude:: out.csv

|

We plot the centerline and plume boundary using
`matplotlib <https://matplotlib.org/>`_

.. plot::
    :context:
    :include-source:

    import matplotlib.pyplot as plt
    import pandas as pd

    df = pd.read_csv("out.csv")

    x = df.x.values
    center = -df.z.values
    upper = -df.z.values + df.radius.values
    lower = -df.z.values - df.radius.values

    plt.plot(x, center, color='k', linewidth=2, label='Centerline')
    plt.plot(x, upper, color='k', linewidth=.5)
    plt.plot(x, lower, color='k', linewidth=.5)
    plt.fill_between(x, lower, upper, color="#e0e0e0",
                     label='Plume extent')

    plt.xlabel('Distance from pipe outlet (m)')
    plt.ylabel('Depth below surface (m)')
    plt.gca().set_aspect('equal')
    plt.legend()
    plt.tight_layout()

|
