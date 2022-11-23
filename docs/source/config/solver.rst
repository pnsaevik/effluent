=========================
Solver parameters
=========================

Solver parameters are technical parameters that modify the solver procedure.
They are passed directly to
`scipy.integrate.solve_ivp <https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html>`_.
All of the parameters have default values which are used unless specified
otherwise.

|

.. confval:: solver.method

    :type: string
    :default: "RK45"

    Integration method to use. Alternatives are "RK45", "RK23, "DOP853",
    "RADAU", "BDF" and "LSODA".

|

.. confval:: solver.rtol

    :type: number
    :default: 1e-3

    Relative tolerance. The solver keeps the local error estimates less than
    atol + rtol * abs(y).

|

.. confval:: solver.atol

    :type: number
    :default: 1e-6

    Absolute tolerance. The solver keeps the local error estimates less than
    atol + rtol * abs(y).

|

.. confval:: solver.first_step

    :type: number
    :default: 0

    Initial step size. Default is 0 which means that the algorithm should choose.

|

.. confval:: solver.max_step

    :type: number
    :default: 0

    Maximum allowed step size. Default is 0 which means that the step size is
    not bounded and determined solely by the solver.
