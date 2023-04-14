=======================================
Buoyant jet in stratified ambient
=======================================

In this example, a horizontal freshwater jet is released into a stratified
water column. The ``config.toml`` file looks like this:

.. literalinclude:: config.toml
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

    # Define unit tangent vector
    xz = np.stack([df.x.values, -df.z.values])
    t = np.roll(xz, shift=-1) - np.roll(xz, shift=1)
    t[:, 0] = xz[:, 1] - xz[:, 0]
    t[:, -1] = xz[:, -1] - xz[:, -2]
    t = t / np.sqrt(np.sum(t**2, axis=0))

    # Define plume boundaries
    x1 = df.x.values - df.radius.values * t[1, :]
    x2 = df.x.values + df.radius.values * t[1, :]
    z1 = -df.z.values + df.radius.values * t[0, :]
    z2 = -df.z.values - df.radius.values * t[0, :]
    x = np.concatenate([x1, np.flip(x2)])
    z = np.concatenate([z1, np.flip(z2)])

    # Generate figure
    plt.plot(xz[0], xz[1], color='k', linewidth=2, label='Centerline')
    plt.fill(x, z, edgecolor='k', linewidth=.5, facecolor="#e0e0e0",
             label='Plume extent')
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
