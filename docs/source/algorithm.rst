===================
Algorithm
===================

Our derivation closely follows |lee2003|_.
The jet discharge is divided into small computational units, or fluid elements.
Each element is a thin cross-sectional slice of the jet which moves and expands
with the flow. The boundaries of the element are constructed so that

1.  The mass and momentum fluxes through the cross-sectional faces of the
    element are zero (the element moves with the flow).

2.  The flux of tracer mass through the lateral boundary of the
    element is zero (the element expands with the jet).

We also assume that there is no net torque acting on the element, and that the
ambient flow is homogeneous at the boundaries of the element. The element will
therefore never rotate.


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
results, we use this formulation in the following derivation. The
end results can be converted back to a gaussian distribution if required.


Jet expansion rate
==================

A turbulent jet expands as it moves through the ambient fluid due to
the entrainment of surrounding water masses at the edges of the jet.
Experiments show that the jet quickly develop a
gaussian profile, which expands laterally at a rate proportional to its speed.
Using the top hat profile, we can express this relation as

.. math ::

    \frac{dR}{dt} = \beta_t \Delta u_t + \beta_n \Delta u_n,

where :math:`\beta_t` is determined by :confval:`model.beta_t`,
:math:`\beta_n` is determined by :confval:`model.beta_n`,
:math:`R` is the jet radius, :math:`t` is time, :math:`\Delta u_t`
is the difference between jet velocity and ambient velocity in the tangential
(along-jet) direction, and :math:`\Delta u_n` is the velocity difference in
the normal (across-jet) direction.

Conservation of mass
====================

Mass increase inside the computational element due to entrainment of ambient
water masses can be expressed as

.. math ::

    \frac{d}{dt}(\rho V) = \rho_a \frac{dV}{dt},

where
:math:`\rho` is the jet density and
:math:`\rho_a` is the ambient water density.


Conservation of momentum
=========================

Change in horizontal momentum due to entrainment of ambient water masses can be
expressed as

.. math ::

    \frac{d}{dt}(\rho u V) =&\, \rho_a u_a \frac{dV}{dt}

    \frac{d}{dt}(\rho v V) =&\, \rho_a v_a \frac{dV}{dt}

where :math:`u` is the horizontal velocity in the direction of the pipe and
:math:`v` is the horizontal transverse velocity with positive direction to the
right of :math:`u`. The subscript :math:`a` denotes ambient quantities.

We assume that the ambient vertical velocity is zero. Vertical momentum change
due to gravity is expressed as

.. math ::

    \frac{d}{dt}(\rho w V) = V K (\rho - \rho_a) g

where :math:`w` is the vertical velocity with positive direction downwards and
:math:`g` is the acceleration of gravity. :math:`K` is the added mass
coefficient, which is a scaling term that reduces the
effect of gravity. The term is required since vertical acceleration of the
plume also stirs up motion of water outside the plume, slowing down the
acceleration. The term depends on the inclination angle of the jet,

.. math ::

    K = k_n \frac{u^2 + v^2}{u^2 + v^2 + w^2} + k_t \frac{w^2}{u^2 + v^2 + w^2},

where :math:`k_n` is determined by :confval:`model.mass_n`
and :math:`k_t` is determined by :confval:`model.mass_t`.

Conservation of volume
=======================

By continuity, the thickness :math:`s` of the computational element is
proportional to the faceward velocity :math:`u`. The volume :math:`V` of the
element can therefore be expressed as

.. math ::

    V = \frac{s_0}{u_0} u \pi R^2,

where the subscript :math:`0` denote initial quantities. The rate of volume
change can be expressed as

.. math ::

    \frac{1}{V} \frac{dV}{dt} = \frac{1}{u} \frac{du}{dt} + \frac{2}{R} \frac{dR}{dt}.


Solving the equations
======================

We choose the following as our primary variables:

==============  =============================================================
Variable        Description
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

The differential equations are reformulated in terms of the primary variables,
and the remaining variables are computed from the primary variables. Using the
vector forms

.. math ::
    \mathbf{x} = x\mathbf{i} + y\mathbf{j} + z\mathbf{k}

and

.. math ::
    \mathbf{u} = u\mathbf{i} + v\mathbf{j} + w\mathbf{k},

we can write the primary equations as:

Displacement
---------------

.. math ::

    \tag{1} \frac{d\mathbf{x}}{dt} = \mathbf{u}

Conservation of momentum:
--------------------------

.. math ::

    \tag{2} \frac{d\mathbf{u}}{dt} = \frac{1}{V} \frac{dV}{dt}  \frac{\rho_a}{\rho} (\mathbf{u}_a - \mathbf{u}) + \frac{1}{\rho} K (\rho - \rho_a) \mathbf{g}

Conservation of mass
------------------------

.. math ::

    \tag{3} \frac{d\rho}{dt} = \frac{1}{V} \frac{dV}{dt} (\rho_a - \rho)

Jet expansion rate
---------------------

.. math ::

    \tag{4} \frac{dR}{dt} = \beta_t \Delta u_t + \beta_n \Delta u_n


The equations are solved using
`scipy.integrate.solve_ivp <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html>`_,
with configurable :doc:`solver parameters </config/solver>`.

Bibliography
===================

.. |lee2003| replace:: Lee and Chu (2003)
.. _lee2003: https://doi.org/10.1007/978-1-4615-0407-8

Lee, Joseph H. W., and Chu, Vincent H. (2003). *Turbulent Jets and Plumes*.
Springer New York, NY.
`doi:10.1007/978-1-4615-0407-8 <https://doi.org/10.1007/978-1-4615-0407-8>`_.
