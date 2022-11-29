===============================
Model parameters
===============================

Model parameters describe different aspects of the
plume model. All of these parameters have default values which are used
unless specified otherwise.

|

.. confval:: model.beta_n

    :type: number
    :default: 0.34

    Entrainment rate coefficient in the normal (across-jet) direction. The
    coefficient is specified relative to the top-hat radius of the jet.

    The default value is from |lee2003|_.

|

.. confval:: model.beta_t

    :type: number
    :default: 0.17

    Entrainment rate coefficient in the tangential (along-jet) direction. The
    coefficient is specified relative to the top-hat radius of the jet.

    The default value is from |lee2003|_.

|

.. confval:: model.mass_n

    :type: number
    :default: 1.0

    Added mass coefficient when gravity acts orthogonal to the jet direction.

    The default value is from |lee2003|_.

|

.. confval:: model.mass_t

    :type: number
    :default: 0.18

    Added mass coefficient when gravity acts tangential to the jet direction.

    The default value is from |lee2003|_.


Bibliography
===================

.. |lee2003| replace:: Lee and Chu (2003)
.. _lee2003: https://doi.org/10.1007/978-1-4615-0407-8

Lee, Joseph H. W., and Chu, Vincent H. (2003). *Turbulent Jets and Plumes*.
Springer New York, NY.
`doi:10.1007/978-1-4615-0407-8 <https://doi.org/10.1007/978-1-4615-0407-8>`_.
