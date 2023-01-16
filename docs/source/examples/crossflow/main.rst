=======================================
Horizontal jet in crossflow
=======================================

A horizontal jet is released into a current, which is in the direction normal
to the pipe outlet. The ambient current velocity is zero initially, but
increases to 0.1 m/s after one hour.
The ``config.toml`` file looks like this:

.. literalinclude:: config.toml
    :language: toml

|

The contents of the output file ``out.csv`` is:

.. literalinclude:: out.csv

|

To visualize the effect of the ambient crossflow, we plot a horizontal
cross-section of the plume at two points in time: Both initially (when the
ambient is stagnant) and after one hour (when there is an ambient current).

.. plot::
    :context: reset
    :include-source:

    import matplotlib.pyplot as plt
    import pandas as pd

    df_all = pd.read_csv("out.csv").groupby('release_time')
    dates = ['1970-01-01 01:00:00', '1970-01-01 02:00:00']
    fcolors = ["#e0707070", "#7070e070"]
    ecolors = ["#a00000", "#0000a0"]

    for date, fcolor, ecolor in zip(dates, fcolors, ecolors):
        df = df_all.get_group(date)

        # Compute unit tangent vector
        velocity = np.sqrt(df.u.values**2 + df.v.values**2)
        tx = df.u.values / velocity
        ty = df.v.values / velocity

        # Compute trajectory of centerline and plume edges
        r = df.radius.values
        x = df.x.values
        y = df.y.values
        x1 = x + r * ty
        y1 = y - r * tx
        x2 = x - r * ty
        y2 = y + r * tx

        plt.plot(
            x, y,
            color=ecolor,
            linewidth=2,
            label=f'{date[-8:-3]} Centerline',
        )
        plt.fill(
            list(x1) + list(reversed(x2)),
            list(y1) + list(reversed(y2)),
            color=fcolor,
            label=f'{date[-8:-3]} Extent',
        )

    plt.xlabel('X distance from pipe outlet (m)')
    plt.ylabel('Y distance from pipe outlet (m)')
    plt.gca().set_aspect('equal')
    plt.legend()
    plt.tight_layout()
|
