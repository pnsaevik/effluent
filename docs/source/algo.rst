===================
Algorithm
===================

Cross-sectional profile
=======================

The jet discharge is divided into small computational units, or fluid elements.
Each element is a thin cross-sectional slice of the jet which moves and expands
with the flow. The boundaries of the element are constructed so that

1.  The mass and momentum fluxes through the cross-sectional faces of the
    element are zero (the element moves with the flow).

2.  The flux of tracer mass through the lateral boundary of the
    element is zero (the element expands with the jet).

In addition, we employ the *top hat profile* simplification as described
by `[1]`_. This means we
assume constant tracer concentration inside the computational element, and
zero outside it, instead of the more realistic gaussian distribution.
It turns out that both the gaussian and top-hat cases are equivalent if

-   The equations are expressed in terms of averages over the computational
    element

-   We require :math:`R = \sigma \sqrt{2}`, where :math:`R` is the radius of
    the element and :math:`\sigma` is the standard deviation of the
    gaussian distribution.

Since the top-hat formulation is easier to work with and yields equivalent
results, we use this formulation extensively in the following derivation. The
end results can be converted back to a gaussian distribution if required.

Jet expansion rate
==================

A turbulent jet expands as it moves through the ambient fluid, due to
the entrainment of surrounding water masses at the edges of the jet.
Experiments show that the concentration and momentum quickly develop a
gaussian profile, which expands laterally at a rate proportional to its speed.
Using the top hat profile, we can express this insight as

.. math ::

    \frac{dR}{dt} = \beta v,

where :math:`R` is the jet radius, :math:`t` is time and :math:`v` is the jet
speed relative to the ambient fluid. The constant :math:`\beta` is determined
by experiments to be :math:`\beta = 0.16`.

Bibliography
===================

.. _[1]:

[1]  Lee, Joseph H. W., and Chu, Vincent H. (2003). *Turbulent Jets and Plumes*.
Boston, MA: Springer US.
`doi:10.1007/978-1-4615-0407-8 <https://doi.org/10.1007/978-1-4615-0407-8>`_.
