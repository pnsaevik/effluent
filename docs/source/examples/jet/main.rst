=======================================
Pure horizontal jet
=======================================

.. plot::
    :context:

    plt.clf()
    plt.close('all')

|

A pure horizontal jet is the simplest test case, which is also investigated in
numerous experimental studies. The ``config.toml`` file looks like this:

.. literalinclude:: config.toml
    :language: toml

|

We run the example from the command line as

.. code-block::

    python -m effluent config.toml

|

or from within python as

.. code-block:: python

    import effluent
    effluent.run("config.toml")  # Alternatively, just pass a dict

|

This produces an output file named ``out.csv``:

.. literalinclude:: out.csv

|

We load the output data using `pandas <https://pandas.pydata.org/>`_

.. plot::
    :context:
    :include-source:

    import pandas as pd
    df = pd.read_csv("out.csv")

|

and plot the centerline and plume boundary using `matplotlib <https://matplotlib.org/>`_

.. plot::
    :context:
    :include-source:

    import matplotlib.pyplot as plt

    x = df.x.values
    center = df.z.values
    upper = df.z.values + df.radius.values
    lower = df.z.values - df.radius.values

    plt.plot(x, center, color='k', linewidth=2, label='Centerline')
    plt.plot(x, upper, color='k', linewidth=.5)
    plt.plot(x, lower, color='k', linewidth=.5)
    plt.fill_between(x, lower, upper, color="#e0e0e0",
                     label='Plume extent')

    plt.xlabel('Distance from pipe outlet (m)')
    plt.ylabel('Depth below surface (m)')
    plt.gcf().set_size_inches(8, 3)
    plt.gca().set_aspect('equal')
    plt.legend()
    plt.tight_layout()

|

Tracer concentration is proportional to plume speed and cross-sectional
area. Here is one way of visualizing the concentration evolution:

.. plot::
    :context:
    :include-source:

    # Extract data
    x = df.x.values
    z = df.z.values
    r = df.radius.values
    element_volume = df.u.values * (r ** 2)
    dilution = element_volume / element_volume[0]

    # Plot plume dilution factor
    xx, zz = np.meshgrid(
        np.linspace(x[0], x[-1], 100),
        np.linspace(z[-1] - r[-1], z[-1] + r[-1], 100),
    )
    rr = np.interp(xx, x, r)
    dilu = np.interp(xx, x, dilution)
    dilu[(z[0] - rr > zz) | (zz > z[0] + rr)] = np.nan
    plt.clf()
    plt.pcolormesh(xx, zz, dilu, cmap='gray', clim=(1, 7))

    # Set visual plot properties
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Dilution factor')
    plt.xlabel('Distance from pipe outlet (m)')
    plt.ylabel('Depth below surface (m)')
    plt.gca().set_aspect('equal')
    plt.gcf().set_size_inches(8, 3)

|

The plume edges are defined using the top-hat profile, as described in
:doc:`../../algorithm`. If desired, we can visualize the plume as a gaussian
profile instead. In a round gaussian plume, the centerline concentration is
twice as large as the corresponding top-hat mean concentration. The top-hat
radius equals the gaussian radius times the square root of two.

.. plot::
    :context:
    :include-source:

    import numpy as np
    plt.clf()

    # Plot fuzzy plume
    depth = df.z.values[0]
    cm = 2 * np.interp(xx, x, 1 / dilution)
    conc = cm * np.exp(-(zz - depth)**2 / (0.5 * rr**2))
    dilu = 1 / (conc + 1e-7)
    plt.pcolormesh(xx, zz, dilu, cmap='gray', clim=(1, 20))

    # Plot centerline and edges
    rad_gauss = df.radius.values / np.sqrt(2)
    lgauss = df.z.values - rad_gauss
    ugauss = df.z.values + rad_gauss
    plt.plot(x, lgauss, color='r', linewidth=.5, label='Gaussian radius')
    plt.plot(x, ugauss, color='r', linewidth=.5)
    plt.plot(x, upper, color='k', linewidth=.5, label='Top-hat radius')
    plt.plot(x, lower, color='k', linewidth=.5)

    # Set visual plot properties
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Dilution factor')
    plt.xlabel('Distance from pipe outlet (m)')
    plt.ylabel('Depth below surface (m)')
    plt.gca().set_aspect('equal')
    plt.gcf().set_size_inches(8, 3)
    plt.legend()

The assumption of gaussian distribution is invalid close to the plume outlet.
A correct representation requires a model of the *Zone of Flow Establishment*.
This is currently not implemented.
