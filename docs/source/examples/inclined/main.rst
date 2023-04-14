=======================================
Inclined heavy jet
=======================================

An inclined jet of salty water enters an environment of
stagnant freshwater.
The ``config.toml`` file looks like this:

.. literalinclude:: config.toml
    :language: toml

|

The contents of the output file ``out.csv`` is

.. literalinclude:: out.csv

|

We plot the centerline and plume boundary using
`matplotlib <https://matplotlib.org/>`_

.. plot::
    :context: reset
    :include-source:

    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np

    df = pd.read_csv("out.csv")

    # Compute tangent vector
    vel = np.sqrt(df.u.values**2 + df.w.values**2)
    tx = df.u.values / vel
    tz = df.w.values / vel

    # Compute plume boundaries
    x1 = df.x.values - df.radius.values * tz
    x2 = df.x.values + df.radius.values * tz
    z1 = -df.z.values - df.radius.values * tx
    z2 = -df.z.values + df.radius.values * tx
    x = np.concatenate([x1, np.flip(x2)])
    z = np.concatenate([z1, np.flip(z2)])

    # Generate figure
    plt.plot(df.x.values, -df.z.values, color='k', linewidth=2, label='Centerline')
    plt.fill(x, z, edgecolor='k', linewidth=.5, facecolor="#e0e0e0", label='Plume extent')
    plt.xlabel('Distance from pipe outlet (m)')
    plt.ylabel('Depth below surface (m)')
    plt.gca().set_aspect('equal')
    plt.legend(loc='lower left')
    plt.tight_layout()

|
