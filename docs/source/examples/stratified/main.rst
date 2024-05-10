=======================================
Buoyant jet in stratified ambient
=======================================

In this example, a horizontal freshwater jet is released into a stratified
water column. The properties of the pipe outflow and the ambient water
masses are specified in external csv files. The ``config.toml`` file looks
like this:

.. literalinclude:: config.toml
    :language: toml

|

Contents of ``pipe.csv`` is shown below. For simplicity, this example
has constant outflow rates. More lines can be added if the outflow
properties are changing with time.

.. literalinclude:: pipe.csv
    :language: toml

|

Contents of ``ambient.csv`` is shown below. More lines can be added
to specify further stratification levels, or to specify properties
that change with time. Sorting order is not critical, but
we recommend sorting by depth, then by time.

.. literalinclude:: ambient.csv
    :language: toml

|

The contents of the output file ``out.csv`` is

.. literalinclude:: out.csv

Observe that integration stopped before the specified end time was
reached. This is because the plume at this point has slowed down so much that
we have entered the far-field regime, as explained in
:ref:`Algorithm <conservation-of-volume>`.

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
    plt.legend(loc='upper left')
    plt.tight_layout()

As seen in the figure, the plume is lifted upwards by buoyancy forces.
At some point, the plume is diluted so much that its buoyancy is neutral
compared to the ambient water masses. It still continues to rise for some time
due to its momentum, overshooting the depth level of neutral buoyancy.
Eventually, it sinks back into a stable depth level.
