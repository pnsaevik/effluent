=======================================
Pure horizontal jet
=======================================

A pure horizontal jet is the simplest test case, which is also investigated in
numerous experimental studies. The ``config.toml`` file looks like this:

.. literalinclude:: config.toml
    :language: toml

|

We run the example from the command line as

.. code-block::

    python -m effluent config.toml

or from within python as

.. code-block:: python

    import effluent
    effluent.run("config.toml")  # Alternatively, just pass a dict

|

We load the output data using ``pandas``

.. plot::
    :context:
    :include-source:

    import pandas as pd
    df = pd.read_csv("out.csv")

and plot the centerline and radius using ``matplotlib``

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
    plt.gca().set_aspect('equal')
    plt.legend()

The plume edges are defined using the top-hat profile, as described in
:doc:`../../algorithm`. If desired, we can visualize the plume as a gaussian
profile instead:


.. plot::
    :context:
    :include-source:

    import numpy as np

    rad_gauss = df.radius.values / np.sqrt(2)
    low_gauss = df.z.values - rad_gauss
    upp_gauss = df.z.values + rad_gauss

    # Plot fuzzy plume
    xx, zz = np.meshgrid(x, np.linspace(lower[-1], upper[-1], 10))
    dist_from_center = zz - df.z.values[0]
    plume = np.exp(-(dist_from_center/rad_gauss)**2)
    plt.cla()
    plt.pcolormesh(xx, zz, plume, cmap="gray_r", shading='gouraud')

    # Plot centerline and edges
    plt.plot(x, low_gauss, color='r', linewidth=.5, label='Gaussian radius')
    plt.plot(x, upp_gauss, color='r', linewidth=.5)
    plt.plot(x, upper, color='k', linewidth=.5, label='Top-hat radius')
    plt.plot(x, lower, color='k', linewidth=.5)

    plt.xlabel('Distance from pipe outlet (m)')
    plt.ylabel('Depth below surface (m)')
    plt.legend()

The assumption of gaussian distribution is invalid close to the plume outlet.
A correct representation requires a model of the *Zone of Flow Establishment*.
This is currently not implemented.
