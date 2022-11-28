import numpy as np
from scipy.integrate import solve_ivp
import xarray as xr
import logging


logger = logging.getLogger(__name__)


class Model:
    def __init__(self, fname_or_dict):
        from effluent.io import load_config, Pipe, Ambient, Output

        self.conf = load_config(fname_or_dict)
        self.pipe = Pipe.from_config(self.conf['pipe'])
        self.ambient = Ambient.from_config(self.conf['ambient'])
        self.output = Output.from_config(self.conf['output'])
        self.solver = Solver.from_config(self.conf['solver'])

    def run(self):
        frequency = self.conf['timestepper']['frequency']
        stop = self.conf['timestepper']['stop']
        times = np.arange(0, stop + frequency / 2, frequency)

        with self.output as output:
            for time in times:
                pipe = self.pipe.select(time)
                ambient = self.ambient.select(time)
                result = self.solver.solve(pipe, ambient)
                output.write(time, result)


class Solver:
    def __init__(self, steps):
        self.steps = steps

    @staticmethod
    def from_config(conf):
        return Solver(**conf)

    def solve(self, pipe, ambient):
        ivp = InitialValueProblem(self.steps, pipe, ambient)
        return ivp.solve()


class InitialValueProblem:
    def __init__(self, steps, pipe, ambient):
        self.steps = steps
        self.pipe = pipe
        self.ambient = ambient
        self._zmin = self.ambient.depth[0].values
        self._zmax = self.ambient.depth[1].values
        self.varnames = ['x', 'y', 'z', 'u', 'v', 'w', 'density', 'radius']
        self.method = 'RK45'

    def initial_conditions(self):
        x0 = 0
        y0 = 0
        z0 = self.pipe.depth.values.item()
        u0 = self.pipe.u.values.item()
        v0 = 0
        w0 = self.pipe.w.values.item()
        d0 = self.pipe.dens.values.item()
        r0 = 0.5 * self.pipe.diam.values.item()
        return np.array([x0, y0, z0, u0, v0, w0, d0, r0], 'f8')

    def solve(self):
        result = solve_ivp(
            fun=self.odefunc,
            t_span=self.steps[[0, -1]],
            y0=self.initial_conditions(),
            method=self.method,
            t_eval=self.steps,
            vectorized=True,
        )

        # noinspection PyUnresolvedReferences
        res_t, res_y = result.t, result.y

        # Organize result
        data_vars = {v: xr.Variable('t', res_y[i]) for i, v in enumerate(self.varnames)}
        return xr.Dataset(data_vars=data_vars, coords=dict(t=res_t))

    def odefunc(self, t, y):
        """
        ODE function to be solved by scipy methods

        The order of the variables is (x, y, z, u, v, w, rho, R)

        :param t: Vectorized time parameter of shape (n_times, )
        :param y: Input vector of shape (n_vars, n_times)
        """

        # Rename input variables
        # noinspection PyUnusedLocal
        t = t
        y_in = y
        x, y, z, u, v, w, rho, R = y_in

        # Define coefficients
        beta_t = 0.16   # Entrainment coefficient, co-flow
        beta_n = 0.4    # Entrainment coefficient, cross-flow
        K_t = 0.85      # Added mass coefficient, tangential gravity pull (= 1 / [1 + k_t])
        K_n = 0.5       # Added mass coefficient, normal gravity pull (= 1 / [1 + k_n])

        # Extract ambient velocity and density
        clipped_depth = np.clip(z, self._zmin, self._zmax)
        ambient = self.ambient.interp(depth=clipped_depth)
        u_a = ambient.u.values
        v_a = ambient.v.values
        rho_a = ambient.dens.values

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

        # Displacement
        ddt_x = u
        ddt_y = v
        ddt_z = w

        # Jet expansion rate (entrainment rate)
        ddt_R = beta_t * delta_u_t + beta_n * delta_u_n

        # Conservation of volume
        ddt_log_R2 = 2 * ddt_R / R
        rho_ratio = rho_a / rho
        ddt_log_V = ddt_log_R2 * u / (u + rho_ratio * delta_u)

        # Conservation of mass
        ddt_rho = ddt_log_V * (rho_a - rho)

        # Conservation of momentum
        prefix = -ddt_log_V * rho_ratio
        ddt_u = prefix * delta_u
        ddt_v = prefix * delta_v
        ddt_w = prefix * w + K * (1 - rho_ratio) * 9.81

        ddt_y = np.stack([ddt_x, ddt_y, ddt_z, ddt_u, ddt_v, ddt_w, ddt_rho, ddt_R])
        return ddt_y
