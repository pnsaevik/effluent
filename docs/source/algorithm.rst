===================
Algorithm
===================

Our derivation closely follows Chapter 9 of |lee2003|_.
The jet discharge is divided into small computational units called fluid elements.
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
It turns out that the conservation equations for the gaussian and top-hat
profiles are equivalent if

-   The equations are expressed in terms of averages over the computational
    element

-   We require :math:`R = \sigma \sqrt{2}`, where :math:`R` is the radius of
    the element and :math:`\sigma` is the standard deviation of the
    gaussian distribution.

Since the top-hat formulation is easier to work with and yields equivalent
results, we use this formulation in the following derivation. The
end results can be converted back to a gaussian distribution if required.

Buoyant jets eventually develop a non-gaussian, kidney shaped cross sectional
profile with a double peak of concentration maxima. In this case, we identify
the plume boundary with the turbulent interface where the intermittency factor
is 50 %.


Jet expansion rate
==================

A turbulent jet expands as it moves through the ambient fluid due to
the entrainment of surrounding water masses at the edges of the jet.
Entrainment is driven by the velocity difference between the jet and the
ambient fluid, which can be decomposed into a tangential and normal component.
The tangential component of the velocity difference is responsible for
*shear entrainment*, while the normal component is responsible for
*vortex entrainment*. To combine the two processes, we employ the
following entrainment hypothesis,

.. math ::
    :label: entrainment

    \frac{dR}{dt} = \beta_t \Delta u_t + \beta_n \Delta u_n,

where :math:`\beta_t` is determined by :confval:`model.beta_t`,
:math:`\beta_n` is determined by :confval:`model.beta_n`,
:math:`R` is the jet radius, :math:`t` is time, :math:`\Delta u_t`
is the difference between jet velocity and ambient velocity in the tangential
(along-jet) direction,

.. math ::
    :label: delta_t

    \Delta u_t = \left| \sqrt{u^2 + v^2 + w^2} - \frac{uu_a+vv_a}{\sqrt{u^2 + v^2 + w^2}} \right|,

and :math:`\Delta u_n` is the velocity difference in the normal (across-jet)
direction,

.. math ::
    :label: delta_n

    \Delta u_n = \sqrt{(u - u_a)^2 + (v - v_a)^2 + w^2 - \Delta u_t^2}.

Conservation of mass
====================

Mass increase inside the computational element due to entrainment of ambient
water masses can be expressed as

.. math ::
    :label: masscons

    \frac{d}{dt}(\rho V) = \rho_a \frac{dV}{dt},

where
:math:`\rho` is the jet density and
:math:`\rho_a` is the ambient water density.


Conservation of momentum
=========================

Change in horizontal momentum due to entrainment of ambient water masses can be
expressed as

.. math ::
    :label: momcons_u

    \frac{d}{dt}(\rho u V) = \rho_a u_a \frac{dV}{dt},

.. math ::
    :label: momcons_v

    \frac{d}{dt}(\rho v V) = \rho_a v_a \frac{dV}{dt},

where :math:`u` is the horizontal velocity in the direction of the pipe and
:math:`v` is the horizontal transverse velocity with positive direction to the
right of :math:`u`. The subscript :math:`a` denotes ambient quantities.

We assume that the ambient vertical velocity is zero. Vertical momentum change
due to gravity is expressed as

.. math ::
    :label: momcons_w

    \frac{d}{dt}(\rho w V) = V K (\rho - \rho_a) g,

where :math:`w` is the vertical velocity with positive direction downwards and
:math:`g` is the acceleration of gravity. :math:`K` is the added mass
coefficient, which is a scaling term that reduces the
effect of gravity. The term is required since vertical acceleration of the
plume also stirs up motion of water outside the plume, slowing down the
acceleration. The term depends on the inclination angle of the jet,

.. math ::
    :label: addmass

    K = \frac{1}{u^2 + v^2 + w^2}\left( \frac{u^2 + v^2}{1 + k_n} + \frac{w^2}{1 + k_t} \right),

where :math:`k_n` is determined by :confval:`model.mass_n`
and :math:`k_t` is determined by :confval:`model.mass_t`.

Conservation of volume
=======================

By continuity, the thickness of the computational element is
proportional to the faceward velocity :math:`u`. The volume :math:`V` of the
element can therefore be expressed as

.. math ::
    :label: voldef

    V = \frac{s_0}{u_0} u \pi R^2,

where :math:`s_0` is the initial thickness and :math:`u_0` the initial
velocity.

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

The differential equations are reformulated in terms of the primary variables.
Remaining variables are computed from the primary variables. Using the
vector forms

.. math ::
    :label: xvec

    \mathbf{x} = x\mathbf{i} + y\mathbf{j} + z\mathbf{k}

and

.. math ::
    :label: uvec

    \mathbf{u} = u\mathbf{i} + v\mathbf{j} + w\mathbf{k},

we can write the primary equations as:

Displacement
---------------

.. math ::
    :label: sol_displ

    \frac{d\mathbf{x}}{dt} = \mathbf{u}

Conservation of momentum:
--------------------------

.. math ::
    :label: sol_mom

    \frac{d\mathbf{u}}{dt} = \frac{1}{V} \frac{dV}{dt}  \frac{\rho_a}{\rho} (\mathbf{u}_a - \mathbf{u}) + \frac{1}{\rho} K (\rho - \rho_a) \mathbf{g}

Conservation of mass
------------------------

.. math ::
    :label: sol_mass

    \frac{d\rho}{dt} = \frac{1}{V} \frac{dV}{dt} (\rho_a - \rho)

Jet expansion rate
---------------------

.. math ::
    :label: sol_jet

    \frac{dR}{dt} = \beta_t \Delta u_t + \beta_n \Delta u_n

|

In addition we utilize the following expression for the rate of volume change,
which is derived from :eq:`masscons`, :eq:`momcons_u` and :eq:`voldef`:

.. math ::
    :label: sol_voldef

    \frac{1}{V}\frac{dV}{dt} = \frac{1}{R}\frac{dR}{dt}\frac{2 \rho u}{\rho u + \rho_a (u - u_a)}


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
