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

We also assume that the element is always oriented in the general direction
of the jet.

Nomenclature
==================

The table below contains symbols that are used on the page. Fluid properties
describe average values within the computational element, unless otherwise
specified.

===================  =================================================================
Variable             Description
===================  =================================================================
:math:`\beta_n`      Normal entrainment rate coefficient, :confval:`model.beta_n`
:math:`\beta_t`      Tangential entrainment rate coefficient, :confval:`model.beta_t`
:math:`\rho`         Mass density
:math:`\rho_a`       Ambient mass density
:math:`K`            Added mass coefficient
:math:`k_n`          Normal added mass coefficient, :confval:`model.mass_n`
:math:`k_t`          Tangential added mass coefficient, :confval:`model.mass_t`
:math:`R`            Radius of the computational element
:math:`\mathbf u`    Three-dimensional vector-valued velocity
:math:`u`            Velocity in the :math:`x` direction
:math:`\mathbf u_a`  Three-dimensional vector-valued ambient velocity
:math:`u_a`          Ambient current velocity in the :math:`x` direction
:math:`\Delta u_n`   Velocity difference, normal (across-jet) direction
:math:`\Delta u_t`   Velocity difference, tangential (along-jet) direction
:math:`V`            Volume of computational element
:math:`v`            Velocity in the :math:`y` direction
:math:`v_a`          Ambient current velocity in the :math:`y` direction
:math:`w`            Velocity in the :math:`z` direction
:math:`\mathbf x`    Three-dimensional vector-valued centerline position of jet
:math:`x`            Horizontal distance from outlet, in the direction parallel to
                     the pipe
:math:`y`            Horizontal distance from outlet, in the direction directly to
                     the right
:math:`z`            Depth below sea surface
===================  =================================================================


Coordinate system
=======================

We choose a coordinate system that is aligned with the outflow pipe. The
:math:`x` variable is the horizontal distance from the outlet in the direction
*parallel* to the pipe. The :math:`y` variable is the horizontal distance from
the outlet in the direction directly to the *right* of the pipe. Finally, the
:math:`z` variable is the depth below sea surface. The jet velocities in each
of these directions are :math:`u`, :math:`v` and :math:`w`, respectively.

The jet centerline position is written in vector form as

.. math ::
    :label: xvec

    \mathbf{x} = x\mathbf{i} + y\mathbf{j} + z\mathbf{k}

and the jet velocity is written in vector form as

.. math ::
    :label: uvec

    \mathbf{u} = u\mathbf{i} + v\mathbf{j} + w\mathbf{k}.


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
end results can be converted back to a gaussian distribution if desired.

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

In the above expressions, :math:`u_a` and :math:`v_a` are the ambient ocean
current velocities in the :math:`u` and :math:`v` directions, respectively.

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
proportional to the faceward velocity. The volume :math:`V` of the
element is therefore given by

.. math ::
    :label: voldef_prim

    \frac{V}{V_0} = \frac{R^2}{R_0^2} \frac{\sqrt{u^2 + v^2 + w^2}}{u_0} ,

where :math:`V_0`, :math:`R_0` and :math:`u_0` are the initial
volume, radius and velocity, respectively. The expression above also represents
the dilution rato of any substance transported by the jet.

By differentiation, we obtain the equivalent equation

.. math ::
    :label: voldef

    \frac{1}{V}\frac{dV}{dt} = 2 \frac{1}{R}\frac{dR}{dt} + \frac{u\frac{du}{dt} + v\frac{dv}{dt} + w\frac{dw}{dt}}{u^2 + v^2 + w^2} .

The volume of the element should always increase. If the expression above is
negative, the jet has slowed down so much that our entrainment rate assumption
is no longer valid. In other words, we have then reached the far-field
dispersion regime. Integration is stopped whenever this happens.

Solving the equations
======================

We choose :math:`\mathbf x`, :math:`\mathbf u`, :math:`\rho` and :math:`R` as
our primary variables. Rewriting the primary equations in terms of these
variables, we obtain:

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
which is derived from :eq:`voldef` and :eq:`sol_mom`:

.. math ::
    :label: sol_voldef

    \frac{1}{V}\frac{dV}{dt}=\frac{\frac{2\rho}{R}\frac{dR}{dt}+K\left(\rho-\rho_{a}\right)\frac{gw}{u^{2}+v^{2}+w^{2}}}{\rho_{a}\left(1-\frac{u_{a}u+v_{a}v}{u^{2}+v^{2}+w^{2}}\right)+\rho}


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
