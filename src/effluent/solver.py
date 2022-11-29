import numpy as np
import xarray as xr
from scipy.integrate import solve_ivp


class Solver:
    def __init__(self):
        self.varnames = ['x', 'y', 'z', 'u', 'v', 'w', 'density', 'radius']

        # Model parameters
        self.beta_n = 0.34
        self.beta_t = 0.17
        self.mass_n = 1.0
        self.mass_t = 0.18

        # Solver parameters
        self.method = "RK45"
        self.rtol = 1e-3
        self.atol = 1e-6
        self.first_step = 0
        self.max_step = 0

        # Output parameters
        self.resolution = None
        self.stagnation = None

        self.data = None

    @staticmethod
    def from_config(conf):
        s = Solver()
        option_names = [
            'beta_n', 'beta_t', 'mass_n', 'mass_t', 'method', 'rtol', 'atol',
            'first_step', 'max_step', 'resolution', 'stagnation',
        ]
        for k, v in conf.items():
            if k in option_names:
                setattr(s, k, v)

        return s

    def solve(self):
        steps = np.arange(0, self.stagnation + 0.5 * self.resolution, self.resolution)

        result = solve_ivp(
            fun=self.odefunc,
            t_span=steps[[0, -1]],
            y0=self.initial_conditions(),
            t_eval=steps,
            vectorized=True,
            method=self.method,
            rtol=self.rtol,
            atol=self.atol,
            first_step=self.first_step or None,
            max_step=self.max_step or np.inf,
        )

        # noinspection PyUnresolvedReferences
        res_t, res_y = result.t, result.y

        # Organize result
        data_vars = {v: xr.Variable('t', res_y[i]) for i, v in enumerate(self.varnames)}
        return xr.Dataset(data_vars=data_vars, coords=dict(t=res_t))

    def ambient_data(self, depth):
        ambient = self.data[1]
        zmin = ambient.depth[0].values
        zmax = ambient.depth[-1].values
        clipped_depth = np.clip(depth, zmin, zmax)
        return ambient.interp(depth=clipped_depth)

    def pipe_data(self):
        return self.data[0]

    def initial_conditions(self):
        pipe = self.pipe_data()
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
        beta_t = self.beta_t          # Entrainment coefficient, co-flow
        beta_n = self.beta_n          # Entrainment coefficient, cross-flow
        K_t = 1 / (1 + self.mass_t)   # Added mass coefficient, tangential gravity pull
        K_n = 1 / (1 + self.mass_n)   # Added mass coefficient, normal gravity pull

        # Extract ambient velocity and density
        ambient = self.ambient_data(z)
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
