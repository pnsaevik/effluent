===================
Algorithm
===================

Our derivation closely follows `[1]`_.
The jet discharge is divided into small computational units, or fluid elements.
Each element is a thin cross-sectional slice of the jet which moves and expands
with the flow. The boundaries of the element are constructed so that

1.  The mass and momentum fluxes through the cross-sectional faces of the
    element are zero (the element moves with the flow).

2.  The flux of tracer mass through the lateral boundary of the
    element is zero (the element expands with the jet).


Cross-sectional profile
=======================

We employ the *top hat profile* simplification. This
means we assume constant tracer concentration inside the computational element
and zero outside it, instead of the more realistic gaussian distribution.
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

A turbulent jet expands as it moves through the ambient fluid due to
the entrainment of surrounding water masses at the edges of the jet.
Experiments show that the jet quickly develop a
gaussian profile, which expands laterally at a rate proportional to its speed.
Using the top hat profile, we can express this relation as

.. math ::

    \frac{dR}{dt} = \beta_1 \Delta u_t + \beta_2 \Delta u_n,

where :math:`R` is the jet radius, :math:`t` is time, :math:`\Delta u_t`
is the difference between jet velocity and ambient velocity in the tangential
(along-jet) direction, and :math:`\Delta u_n` is the velocity difference in
the normal (across-jet) direction. The constants :math:`\beta_1` and
:math:`\beta_2` are determined by experiments to be
:math:`\beta_1 = 0.16` and :math:`\beta_2 = 0.4`.


Conservation of mass
====================

Mass increase inside the computational element is due to entrainment of ambient
water masses. This can be expressed as

.. math ::

    \frac{d}{dt}(A \rho) = \frac{dA}{dt}\rho_a

where :math:`A = \pi R^2` is the area, :math:`\rho` is the jet density and
:math:`\rho_a` is the ambient water density.


Conservation of momentum
=========================

Entrainment of ambient water masses into the computational element also changes
the momentum. This can be expressed as

.. math ::

    \frac{d}{dt}(A \rho u) =&\, \frac{dA}{dt}\rho_a u_a

    \frac{d}{dt}(A \rho v) =&\, \frac{dA}{dt}\rho_a v_a

    \frac{d}{dt}(A \rho w) =&\, A K (\rho - \rho_a) g

where :math:`g` is the acceleration of gravity, :math:`K` is the added mass
coefficient (explained below), :math:`(u, v, w)` is the vector-valued velocity
and the subscript :math:`a` denotes ambient quantities. In vector notation,

.. math ::

    \frac{d}{dt}(A \rho \mathbf{u}) = \frac{dA}{dt}\rho_a \mathbf{u_a} + A K (\rho - \rho_a) \mathbf{g},

In our chosen coordinate system, :math:`u` is the horizontal velocity
in the direction of the pipe, :math:`v` is the horizontal transverse velocity,
with positive direction to the right of :math:`u`, while :math:`w` is the
vertical velocity, with positive direction downwards.

The added mass coefficient :math:`K` is a scaling term which reduces the
effect of gravity, due to the fact that vertical acceleration of the plume also
stirs up motion of water outside the plume. The term is dependent on the
inclination angle of the jet,

.. math ::

    K = k_t \frac{u^2 + v^2}{u^2 + v^2 + w^2} + k_n \frac{w^2}{u^2 + v^2 + w^2}

with :math:`k_t = 0.5` and  :math:`k_t = 0.85`.


Solving the equations
======================

We choose the following variables as our primary variables. The differential
equations are reformulated in terms of the primary variables, and the remaining
variables are computed from the primary variables.

==============  =============================================================
:math:`x`       Horizontal distance from outlet, in the direction parallel to
                the pipe
:math:`y`       Horizontal distance from outlet, in the direction directly to
                the right
:math:`z`       Depth below sea surface
:math:`u`       Velocity in the :math:`x` direction
:math:`v`       Velocity in the :math:`y` direction
:math:`w`       Velocity in the :math:`z` direction
:math:`\rho`     Mass density
:math:`R`       Radius of the computational element
==============  =============================================================

Reformulated equations below:

**Change of area, by definition**

.. math ::

    \frac{1}{A} \frac{dA}{dt} = \frac{2}{b} \frac{db}{dt}

**Conservation of mass**

.. math ::

    \frac{d\rho}{dt} = \frac{1}{A} \frac{dA}{dt} (\rho_a - \rho)

**Conservation of momentum:**

.. math ::

    \frac{d\mathbf{u}}{dt} = \frac{1}{A} \frac{dA}{dt}  \frac{\rho_a}{\rho} (\mathbf{u}_a - \mathbf{u}) + \frac{1}{\rho} K (\rho - \rho_a) \mathbf{g}

Bibliography
===================

.. _[1]:

[1]  Lee, Joseph H. W., and Chu, Vincent H. (2003). *Turbulent Jets and Plumes*.
Boston, MA: Springer US.
`doi:10.1007/978-1-4615-0407-8 <https://doi.org/10.1007/978-1-4615-0407-8>`_.
