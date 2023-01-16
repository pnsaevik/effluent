=======================================
Horizontal jet in coflow
=======================================

.. plot::
    :context:

    plt.close('all')

|

A horizontal jet is released into a current, which is in the same direction
as the pipe outlet. The ambient current velocity is zero initially, but
increases to 0.1 m/s after one hour.
The ``config.toml`` file looks like this:

.. literalinclude:: config.toml
    :language: toml

|

The contents of the output file ``out.csv`` is:

.. literalinclude:: out.csv

|

To visualize the effect of the ambient flow, we plot a horizontal
cross-section of the plume at two points in time: Both initially (when the
ambient is stagnant) and after one hour (when there is an ambient current).

.. plot::
    :context:
    :include-source:

    import matplotlib.pyplot as plt
    import pandas as pd

    df_all = pd.read_csv("out.csv").groupby('release_time')
    dates = ['1970-01-01 01:00:00', '1970-01-01 02:00:00']
    fcolors = ["#e0707070", "#7070e070"]
    ecolors = ["#a00000", "#0000a0"]

    for date, fcolor, ecolor in zip(dates, fcolors, ecolors):
        df = df_all.get_group(date)

        # Compute trajectory of centerline and plume edges
        r = df.radius.values
        x = df.x.values
        y = df.y.values
        y1 = y - r
        y2 = y + r

        plt.plot(
            x, y,
            color=ecolor,
            linewidth=2,
            label=f'{date[-8:-3]} Centerline',
        )
        plt.fill_between(
            x, y1, y2,
            color=fcolor,
            label=f'{date[-8:-3]} Extent',
        )

    plt.xlabel('X distance from pipe outlet (m)')
    plt.ylabel('Y distance from\npipe outlet (m)')
    plt.gca().set_aspect('equal')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.35), ncol=5)
    plt.gcf().set_size_inches(8, 3)
    plt.tight_layout()

|
