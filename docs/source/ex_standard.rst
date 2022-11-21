================================
Standard example
================================

In this example we show how to generate some code. We start with the following:

.. doctest::

    >>> import effluent
    >>> 2 + 3
    5

As we can see, the end result is as expected.

However, the xarray import may be more challenging..

Let us see if it works with matplotlib:

.. plot::

    import effluent
    import matplotlib.pyplot as plt
    import numpy as np
    x = np.random.randn(1000)
    plt.hist( x, 20)
    plt.grid()
    plt.title(f'Version: {effluent.__version__}')
    plt.show()
