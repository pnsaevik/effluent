import numpy as np
from scipy.integrate import solve_ivp
import xarray as xr


class Solver:
    def __init__(self, conf):
        self.conf = conf
        self.stop = conf['stop']
        self.resolution = conf['resolution']
        self.frequency = conf['frequency']
        self.stagnation = conf['stagnation']

    def solve(self, pipe, ambient):

        # Setup equation
        varnames = ['x', 'y', 'z', 'u', 'v', 'w', 'density', 'radius']
        init = np.zeros(8, dtype='f8')
        steps = np.arange(0, self.stagnation, self.resolution)
        fun = lambda t, y: np.zeros_like(y)
        release_time = 0

        # Solve equation
        result = solve_ivp(
            fun=fun,
            t_span=steps[[0, -1]],
            y0=init,
            method='RK45',
            t_eval=steps,
            vectorized=True,
        )

        # noinspection PyUnresolvedReferences
        result_t, result_y = result.t, result.y

        # Organize result
        data_vars = {v: xr.Variable('t', result_y[i]) for i, v in enumerate(varnames)}
        coords = dict(
            release_time=xr.Variable(dims=(), data=release_time),
            t=xr.Variable(dims='t', data=result_t),
        )

        data_vars['x'].attrs['long_name'] = 'centerline x coordinate'
        data_vars['y'].attrs['long_name'] = 'centerline y coordinate'
        data_vars['z'].attrs['long_name'] = 'centerline z coordinate'
        data_vars['u'].attrs['long_name'] = 'velocity in x direction'
        data_vars['v'].attrs['long_name'] = 'velocity in y direction'
        data_vars['w'].attrs['long_name'] = 'velocity in z direction'
        data_vars['density'].attrs['long_name'] = 'mass density of fluid'
        data_vars['radius'].attrs['long_name'] = 'radius of top hat profile'
        coords['t'].attrs['long_name'] = 'time since release'
        coords['release_time'].attrs['long_name'] = 'time of release'

        data_vars['x'].attrs['units'] = 'm'
        data_vars['y'].attrs['units'] = 'm'
        data_vars['z'].attrs['units'] = 'm'
        data_vars['u'].attrs['units'] = 'm/s'
        data_vars['v'].attrs['units'] = 'm/s'
        data_vars['w'].attrs['units'] = 'm/s'
        data_vars['density'].attrs['units'] = 'kg/m^3'
        data_vars['radius'].attrs['units'] = 'm'
        coords['t'].attrs['units'] = 's'
        coords['release_time'].attrs['units'] = 's'

        data_vars['z'].attrs['standard_name'] = 'depth_below_surface'
        data_vars['z'].attrs['positive'] = 'down'

        dset = xr.Dataset(
            data_vars=data_vars,
            coords=coords,
        )

        return dset
