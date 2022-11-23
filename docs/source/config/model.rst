===============================
Model parameters
===============================

Model parameters describe different aspects of the
plume model. All of these parameters have default values which are used
unless specified otherwise.

|

.. confval:: model.beta_n

    :type: number
    :default: 0.4

    Entrainment rate coefficient in the normal (across-jet) direction. The
    coefficient is specified relative to the top-hat radius of the jet.

|

.. confval:: model.beta_t

    :type: number
    :default: 0.16

    Entrainment rate coefficient in the tangential (along-jet) direction. The
    coefficient is specified relative to the top-hat radius of the jet.

|

.. confval:: model.mass_n

    :type: number
    :default: 0.5

    Added mass coefficient when gravity acts orthogonal to the jet direction.

|

.. confval:: model.mass_t

    :type: number
    :default: 0.85

    Added mass coefficient when gravity acts tangential to the jet direction.
