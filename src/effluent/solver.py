"""
The module contains the :class:`Solver <effluent.solver.Solver>` class, which performs
numerical integration of the differential equations described in :doc:`/algorithm`.
"""

import numpy as np
import scipy.integrate
import effluent.io
import xarray as xr
import logging


logger = logging.getLogger(__name__)


class Solver:
    """
    The class contains methods for solving the model equations.

    The constructor contains many adjustable parameters, related to the numerical
    solution procedure, the physical model parametrization and the density of output
    points. A more detailed description of each parameter is given in
    :doc:`/config/solver`, :doc:`/config/model` and :doc:`/config/output`.

    :param beta_n: Entrainment rate coefficient in the normal (across-jet) direction
    :param beta_t: Entrainment rate coefficient in the tangential (along-jet) direction
    :param mass_n: Added mass coefficient in the orthogonal direction
    :param mass_t: Added mass coefficient in the tangential direction
    :param method: Integration method to use
    :param rtol: Relative tolerance
    :param atol: Absolute tolerance
    :param first_step: Initial step size
    :param max_step: Maximum allowed step size
    :param start: Time (in seconds) of first trajectory point
    :param stop: Time (in seconds) of last trajectory point
    :param step: Time (in seconds) between trajectory points
    """

    def __init__(self, beta_n=0.34, beta_t=0.17, mass_n=1.0, mass_t=0.18, method="RK45",
                 rtol=1e-3, atol=1e-6, first_step=0, max_step=0, start=0, stop=60,
                 step=1):
        self.varnames = ['x', 'y', 'z', 'u', 'v', 'w', 'density', 'radius']

        # Model parameters
        self.beta_n = beta_n
        self.beta_t = beta_t
        self.mass_n = mass_n
        self.mass_t = mass_t

        # Solver parameters
        self.method = method
        self.rtol = rtol
        self.atol = atol
        self.first_step = first_step
        self.max_step = max_step

        # Output parameters
        self.start = start
        self.stop = stop
        self.step = step

        self._pipe = None
        self._ambient = None

    def set_init(self, pipe: effluent.io.Pipe, time):
        """
        Set initial conditions

        Set initial conditions of the initial value problem by sampling pipe parameters
        at a specific point in time.

        :param pipe: Properties of the pipe and effluent discharge
        :param time: Sampling time
        """
        self._pipe = pipe.select(time).compute()

    def set_ambient(self, ambient: effluent.io.Ambient, time):
        """
        Set ambient conditions

        Set ambient conditions of the initial value problem by sampling ocean data
        at a specific point in time.

        :param ambient: Properties of the ambient ocean
        :param time: Sampling time
        """
        self._ambient = ambient.select(time).compute()

        # Check that the data is sorted by increasing depth
        depth = self._ambient.depth.values
        if not np.all(depth[:-1] < depth[1:]):
            logger.error("Ambient conditions are not sorted by increasing depth")
            raise ValueError("Ambient conditions are not sorted by increasing depth")

    def volume_change_ratio(self, t, y):
        """
        Compute the time derivative of log(V) according to :eq:`sol_voldef`.

        This function is required to determine when the simulation should be terminated.

        :param t: Seconds since release
        :param y: Tuple of primary variables: ``x``, ``y``, ``z``, ``u``, ``v``, ``w``, ``density``, ``radius``
        :return: Time derivative of log(V)
        """
        _vars = self._compute_vars(t, y)
        return _vars[0]

    def _compute_vars(self, t, y):
        # Rename input variables
        # noinspection PyUnusedLocal
        t = t
        y_in = y
        x, y, z, u, v, w, rho, R = y_in

        # Define coefficients
        beta_t = self.beta_t          # Entrainment coefficient, co-flow
        beta_n = self.beta_n          # Entrainment coefficient, cross-flow
        K_t = 1 / (1 + self.mass_t)   # Added mass coefficient, tangential gravity pull
        K_n = 1 / (1 + self.mass_n)   # Added mass coefficient, normal gravity pull
        g = 9.81                      # Acceleration of gravity

        # Extract ambient velocity and density
        rho_a, u_a, v_a = self._ambient_data(z)

        # Compute added mass coefficient
        squared_horizontal_speed = u*u + v*v
        w2 = w*w
        squared_speed = squared_horizontal_speed + w2
        K = (K_n * squared_horizontal_speed + K_t * w2) / squared_speed

        # Compute flow difference in tangential and normal direction
        delta_u = u - u_a
        delta_v = v - v_a
        squared_excess_speed = delta_u*delta_u + delta_v*delta_v + w2
        speed = np.sqrt(squared_speed)
        delta_u_t = np.abs(speed - (u * u_a + v * v_a) / speed)
        sq_delta_u_n = squared_excess_speed - delta_u_t * delta_u_t
        delta_u_n = np.sqrt(np.maximum(0, sq_delta_u_n))

        # Jet expansion rate (entrainment rate)
        ddt_R = beta_t * delta_u_t + beta_n * delta_u_n

        # Conservation of volume
        ddt_log_R2 = 2 * ddt_R / R
        rho_ratio = rho_a / rho
        gravity_factor = K * (1 - rho_ratio) * g
        nominator = ddt_log_R2 + gravity_factor * w / squared_speed
        denominator = rho_ratio * (1 - (u * u_a + v * v_a) / squared_speed) + 1
        ddt_log_V = nominator / denominator

        return (ddt_log_V, rho_a, rho, rho_ratio, delta_u, delta_v, u, v, w,
                gravity_factor, ddt_R)

    def solve(self) -> xr.Dataset:
        """
        Solve the differential equations described in :doc:`/algorithm`.

        Internally, this function uses ``scipy.integrate.solve_ivp`` to compute the
        solution. Solver properties are set using the class constructor.

        The returned dataset contains the value of each primary variable at the
        pre-determined output points. The variable names are ``x``, ``y``, ``z``, ``u``,
        ``v``, ``w``, ``density``, ``radius``, all indexed by the coordinate ``t``.

        :return: An xarray.Dataset containing the solution
        """
        steps = np.arange(self.start, self.stop + 0.5 * self.step, self.step)

        event = lambda t, y: self.volume_change_ratio(t, y)
        event.terminal = True
        event.direction = -1

        # noinspection PyUnresolvedReferences
        result = scipy.integrate.solve_ivp(
            fun=self.odefunc,
            t_span=steps[[0, -1]],
            y0=self._initial_conditions(),
            t_eval=steps,
            vectorized=True,
            method=self.method,
            rtol=self.rtol,
            atol=self.atol,
            first_step=self.first_step or None,
            max_step=self.max_step or np.inf,
            events=event,
        )

        # noinspection PyUnresolvedReferences
        res_t, res_y, evt_t, evt_y = result.t, result.y, result.t_events, result.y_events

        # Append the end result, if integration was stopped
        if len(evt_t[0]) > 0:
            res_t = np.concatenate([res_t, evt_t[0]])
            res_y = np.concatenate([res_y, evt_y[0].T], axis=1)

        # Organize result
        data_vars = {v: xr.Variable('t', res_y[i]) for i, v in enumerate(self.varnames)}
        return xr.Dataset(data_vars=data_vars, coords=dict(t=res_t))

    def _ambient_data(self, depth):
        amb = self._ambient
        d_a = np.interp(depth, amb.depth.values, amb.dens.values)
        u_a = np.interp(depth, amb.depth.values, amb.u.values)
        v_a = np.interp(depth, amb.depth.values, amb.v.values)

        return np.array([d_a, u_a, v_a])

    def _initial_conditions(self):
        pipe = self._pipe
        x0 = 0
        y0 = 0
        z0 = pipe.depth.values.item()
        u0 = pipe.u.values.item()
        v0 = 0
        w0 = pipe.w.values.item()
        d0 = pipe.dens.values.item()
        r0 = 0.5 * pipe.diam.values.item()
        return np.array([x0, y0, z0, u0, v0, w0, d0, r0], 'f8')

    def odefunc(self, t, y):
        """
        ODE function to be solved by scipy methods

        The order of the variables is (x, y, z, u, v, w, density, radius)

        :param t: Vectorized time parameter of shape (n_times, )
        :param y: Input vector of shape (n_vars, n_times)
        :return: Time derivative of the primary variables (n_vars, n_times)
        """

        _vars = self._compute_vars(t, y)
        ddt_log_V, rho_a, rho, rho_ratio, delta_u, delta_v, u, v, w, gf, ddt_R = _vars

        # Conservation of mass
        ddt_rho = ddt_log_V * (rho_a - rho)

        # Conservation of momentum
        prefix = -ddt_log_V * rho_ratio
        ddt_u = prefix * delta_u
        ddt_v = prefix * delta_v
        ddt_w = prefix * w + gf

        # Displacement
        ddt_x = u
        ddt_y = v
        ddt_z = w

        ddt_y = np.stack([ddt_x, ddt_y, ddt_z, ddt_u, ddt_v, ddt_w, ddt_rho, ddt_R])
        return ddt_y
