"""
The module contains the :class:`Solver <effluent.solver.Solver>` class, which performs
numerical integration of the differential equations described in :doc:`/algorithm`.
"""

import numpy as np
import scipy.integrate
import effluent.io
import xarray as xr
import logging
from collections import namedtuple


logger = logging.getLogger(__name__)
SolverVars = namedtuple('SolverVars', ['x', 'y', 'z', 'u', 'v', 'w', 'density', 'radius', 'tracers'])
OdefuncVars = namedtuple('OdefunVars', ['ddt'])


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
                 rtol=1e-3, atol=1e-6, first_step=0, max_step=0, start=0, stop=60, step=1):
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

        # Internal dataset to hold the initial conditions
        self._pipe = xr.Dataset()

        # Internal dataset to hold the ambient conditions
        self._ambient = xr.Dataset()

        # Internal list of tracked tracer variables
        self._tracers = ('salt', 'temp')

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
        _ = t  # Silence warning about variable not used

        solver_vars = self._unpack(y)
        _, ddt_log_V = self._odefunc(solver_vars)
        return ddt_log_V

    @staticmethod
    def _pack(solver_vars: SolverVars) -> np.ndarray:
        """
        Pack a set of solver variables and tracer variables into an ODE array

        :param solver_vars: A set of solver variables
        :return: Numpy array to be fed into scipy ODE solver
        """
        primary_vars = list(solver_vars)[:-1]
        tracer_vars = list(solver_vars.tracers.values())
        return np.stack(primary_vars + tracer_vars)

    def _unpack(self, ode_array: np.ndarray) -> SolverVars:
        """
        Unpack an ODE array to solver variables and tracer variables

        :param ode_array: A numpy array designed to work with scipy ODE solvers
        :return: A set of solver variables
        """
        num_primary_vars = len(SolverVars._fields) - 1
        primary_vars = ode_array[:num_primary_vars]
        tracer_vars = ode_array[num_primary_vars:]
        tracer_dict = dict(zip(self._tracers, tracer_vars))
        solver_vars = SolverVars(*primary_vars, tracers=tracer_dict)
        return solver_vars

    def _odefunc(self, solver_vars: SolverVars) -> tuple[SolverVars, np.ndarray]:
        """
        ODE function to be used by the solver

        The function implements the formulas described in :doc:`/algorithm`.
        For enhanced readability, the function has annotated arrays (namedtuple
        objects) as input and output. The function is wrapped by self.odefunc
        which packages the input and output as numpy arrays to fit with scipy
        solvers.

        Output from this function includes time derivatives of all solver
        variables as well as auxillary variables used for event detection

        :param solver_vars: Solver variables as described in :doc:`/algorithm`.
        :return: Tuple (odefunc_vars, ddt_log_V), where odefunc_vars are time
            derivatives of solver variables and ddt_log_V is used for event
            detection.
        """

        # Define coefficients
        beta_t = self.beta_t          # Entrainment coefficient, co-flow
        beta_n = self.beta_n          # Entrainment coefficient, cross-flow
        K_t = 1 / (1 + self.mass_t)   # Added mass coefficient, tangential gravity pull
        K_n = 1 / (1 + self.mass_n)   # Added mass coefficient, normal gravity pull
        g = 9.81                      # Acceleration of gravity

        # Extract variables
        z = solver_vars.z
        u = solver_vars.u
        v = solver_vars.v
        w = solver_vars.w

        # Extract ambient velocity and density
        rho_a, u_a, v_a, tracer_ambients = self._ambient_data(z)

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
        ddt_log_R2 = 2 * ddt_R / solver_vars.radius
        rho_ratio = rho_a / solver_vars.density
        gravity_factor = K * (1 - rho_ratio) * g
        nominator = ddt_log_R2 + gravity_factor * w / squared_speed
        denominator = rho_ratio * (1 - (u * u_a + v * v_a) / squared_speed) + 1
        ddt_log_V = nominator / denominator

        # Conservation of mass
        ddt_rho = ddt_log_V * (rho_a - solver_vars.density)

        # Tracers (similar to conservation of mass)
        ddt_tracers = {}
        for tracer_name, tracer_ambient in tracer_ambients.items():
            tracer_value = solver_vars.tracers[tracer_name]
            ddt_tracers[tracer_name] = ddt_log_V * (tracer_ambient - tracer_value)

        # Conservation of momentum
        prefix = -ddt_log_V * rho_ratio
        ddt_u = prefix * delta_u
        ddt_v = prefix * delta_v
        ddt_w = prefix * w + gravity_factor

        # Displacement
        ddt_x = u
        ddt_y = v
        ddt_z = w

        # Assemble output values
        odefunc_vars = SolverVars(
            x=ddt_x,
            y=ddt_y,
            z=ddt_z,
            u=ddt_u,
            v=ddt_v,
            w=ddt_w,
            density=ddt_rho,
            radius=ddt_R,
            tracers=ddt_tracers,
        )

        return odefunc_vars, ddt_log_V

    def solve(self) -> xr.Dataset:
        """
        Solve the differential equations described in :doc:`/algorithm`.

        Internally, this function uses ``scipy.integrate.solve_ivp`` to compute the
        solution. Solver properties are set using the class constructor.

        The returned dataset contains the value of each primary variable at the
        pre-determined output points. The variable names are ``x``, ``y``, ``z``, ``u``,
        ``v``, ``w``, ``density``, ``radius``, all indexed by the coordinate ``t``.

        The returned dataset also contains the secondary variable ``dilution``,
        indexed by the coordinate ``t``.

        :return: An xarray.Dataset containing the solution
        """
        steps = np.arange(self.start, self.stop + 0.5 * self.step, self.step)

        def event(t, y):
            return self.volume_change_ratio(t, y)

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
        data_dict = self._unpack_to_dict(res_y)
        data_vars = {
            k: xr.Variable('t', v) for k, v in data_dict.items()
            # Ignore unavailable tracers
            if k not in self._tracers or k in self._ambient or k in self._pipe
        }
        data_vars['dilution'] = xr.Variable(
            data=self._dilution_factor(res_y),
            dims='t',
        )
        return xr.Dataset(data_vars=data_vars, coords=dict(t=res_t))

    def _unpack_to_dict(self, y: np.ndarray) -> dict:
        """
        Convert an array of scipy ODE style into a dict of named variables

        :param y: Array of scipy ODE style
        :return: Dict of named variables
        """
        names = list(SolverVars._fields)[:-1]
        names += list(self._tracers)
        return dict(zip(names, y))

    def _ambient_data(self, depth):
        amb = self._ambient
        depths = amb.depth.values
        d_a = np.interp(depth, depths, amb.dens.values)
        u_a = np.interp(depth, depths, amb.u.values)
        v_a = np.interp(depth, depths, amb.v.values)

        # Load tracer data, or replace with zero if missing
        t_a = {}
        for tracer_name in self._tracers:
            if tracer_name in amb.variables:
                values = amb[tracer_name].values
                t_a[tracer_name] = np.interp(depth, depths, values)
            else:
                t_a[tracer_name] = np.zeros(np.shape(depth))

        return d_a, u_a, v_a, t_a

    def _initial_conditions(self):
        pipe = self._pipe

        # Load tracer data, or replace with zero if missing
        tracers = {}
        for tracer_name in self._tracers:
            if tracer_name in pipe.variables:
                tracers[tracer_name] = pipe[tracer_name].values.item()
            else:
                tracers[tracer_name] = 0

        init_values = SolverVars(
            x=0,
            y=0,
            z=pipe.depth.values.item(),
            u=pipe.u.values.item(),
            v=0,
            w=pipe.w.values.item(),
            density=pipe.dens.values.item(),
            radius=0.5 * pipe.diam.values.item(),
            tracers=tracers,
        )
        return self._pack(init_values)

    def _dilution_factor(self, y_in):
        """
        Compute dilution factor according to formula :eq:`voldef_prim` in
        :doc:`/algorithm`.

        :param y_in: Solver variables in scipy ODE format
        :return: Dilution factor
        """
        sv0 = self._unpack(self._initial_conditions())
        sv1 = self._unpack(y_in)
        dilution = (
            ((sv1.radius * sv1.radius) / (sv0.radius * sv0.radius)) *
            (np.sqrt(sv1.u * sv1.u + sv1.v * sv1.v + sv1.w * sv1.w) /
             np.sqrt(sv0.u * sv0.u + sv0.v * sv0.v + sv0.w * sv0.w))
        )
        return dilution

    def odefunc(self, t, y):
        """
        ODE function to be solved by scipy methods

        The order of the variables is (x, y, z, u, v, w, density, radius)

        :param t: Vectorized time parameter of shape (n_times, )
        :param y: Input vector of shape (n_vars, n_times)
        :return: Time derivative of the primary variables (n_vars, n_times)
        """
        _ = t  # Silence warning about variable not used

        solver_vars = self._unpack(y)
        odefunc_vars, _ = self._odefunc(solver_vars)
        ddt_y = self._pack(odefunc_vars)
        return ddt_y
