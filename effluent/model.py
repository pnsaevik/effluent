import numpy as np
from scipy.integrate import solve_ivp
import xarray as xr
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


class Model:
    def __init__(self, fname_or_dict):
        self.conf = Config(fname_or_dict)
        pipe, ambient, output, solver = self.conf.init_modules()
        self.pipe = pipe
        self.ambient = ambient
        self.output = output
        self.solver = solver

    def run(self):
        frequency = self.conf.timestepper['frequency']
        stop = self.conf.timestepper['stop']
        times = np.arange(0, stop + frequency / 2, frequency)

        with self.output as output:
            for time in times:
                pipe = self.pipe.select(time)
                ambient = self.ambient.select(time)
                result = self.solver.solve(pipe, ambient)
                output.write(time, result)


class Config:
    def __init__(self, fname_or_dict):
        self.conf = load_config(fname_or_dict)
        self._solver = None
        self._output = None
        self._timestepper = None

    @property
    def solver(self):
        if self._solver is None:
            self._solver = self._generate_solver_conf()
        return self._solver

    @property
    def timestepper(self):
        if self._timestepper is None:
            self._timestepper = self._generate_timestepper_conf()
        return self._timestepper

    @property
    def pipe(self):
        return {}

    @property
    def ambient(self):
        return {}

    @property
    def output(self):
        if self._output is None:
            self._output = self._generate_output_conf()
        return self._output

    def _generate_solver_conf(self):
        out_conf = self.conf['output']
        keys = {'resolution', 'stagnation', 'frequency', 'stop'}
        conf = {k: v for k, v in out_conf.items() if k in keys}
        return conf

    def _generate_output_conf(self):
        out_conf = self.conf['output']
        keys = {'file'}
        conf = {k: v for k, v in out_conf.items() if k in keys}
        return conf

    def _generate_timestepper_conf(self):
        out_conf = self.conf['output']
        keys = {'frequency', 'stop'}
        conf = {k: v for k, v in out_conf.items() if k in keys}
        return conf

    def init_modules(self):
        pipe = Pipe(**self.pipe)
        ambient = Ambient(**self.ambient)
        output = Output.open(**self.output)
        solver = Solver(**self.solver)

        return pipe, ambient, output, solver


def load_config(fname_or_dict):
    if isinstance(fname_or_dict, dict):
        conf = fname_or_dict
    else:
        import yaml
        with open(fname_or_dict, encoding='utf-8') as fp:
            conf = yaml.safe_load(fp)

    return conf


class Pipe:
    def __init__(self):
        pass

    def select(self, time):
        return NotImplemented


class Ambient:
    def __init__(self):
        pass

    def select(self, time):
        return NotImplemented


class Output:
    @staticmethod
    def open(file):
        suffix = Path(file).suffix
        subclasses = {'.csv': OutputCSV, '.nc': OutputNC}
        subclass = subclasses.get(suffix, subclasses['.nc'])
        return subclass(file)


class OutputCSV(Output):
    def __init__(self, file):
        self.file = file
        self._file = None
        self._blank_file = True

    def __enter__(self):
        self._file = open(self.file, 'w', encoding='utf-8', newline='\n')
        self._blank_file = True
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        if self._file is not None:
            self._file.close()
            self._file = None

    def write(self, time, result):
        df = result.to_dataframe()
        df['release_time'] = time
        df = df.reset_index().set_index(['release_time', 't'])

        # Append result to file, write headers only if blank file
        df.to_csv(self._file, line_terminator='\n', header=self._blank_file)
        self._blank_file = False


class OutputNC(Output):
    def __init__(self, file):
        self.file = file

    @staticmethod
    def _append_attributes(result):

        result = result.assign_coords(release_time=time)

        r = result
        r.coords['release_time'].attrs['long_name'] = 'time of release'
        r.coords['release_time'].attrs['units'] = 's'


class Solver:
    def __init__(self, resolution, stagnation, frequency, stop):
        self.stop = stop
        self.resolution = resolution
        self.frequency = frequency
        self.stagnation = stagnation

    def solve(self, pipe, ambient):

        # Setup equation
        varnames = ['x', 'y', 'z', 'u', 'v', 'w', 'density', 'radius']
        init = np.zeros(8, dtype='f8')
        steps = np.arange(0, self.stagnation, self.resolution)
        fun = lambda t, y: np.zeros_like(y)

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
        coords = dict(t=xr.Variable(dims='t', data=result_t))

        data_vars['x'].attrs['long_name'] = 'centerline x coordinate'
        data_vars['y'].attrs['long_name'] = 'centerline y coordinate'
        data_vars['z'].attrs['long_name'] = 'centerline z coordinate'
        data_vars['u'].attrs['long_name'] = 'velocity in x direction'
        data_vars['v'].attrs['long_name'] = 'velocity in y direction'
        data_vars['w'].attrs['long_name'] = 'velocity in z direction'
        data_vars['density'].attrs['long_name'] = 'mass density of fluid'
        data_vars['radius'].attrs['long_name'] = 'radius of top hat profile'
        coords['t'].attrs['long_name'] = 'time since release'

        data_vars['x'].attrs['units'] = 'm'
        data_vars['y'].attrs['units'] = 'm'
        data_vars['z'].attrs['units'] = 'm'
        data_vars['u'].attrs['units'] = 'm/s'
        data_vars['v'].attrs['units'] = 'm/s'
        data_vars['w'].attrs['units'] = 'm/s'
        data_vars['density'].attrs['units'] = 'kg/m^3'
        data_vars['radius'].attrs['units'] = 'm'
        coords['t'].attrs['units'] = 's'

        data_vars['z'].attrs['standard_name'] = 'depth_below_surface'
        data_vars['z'].attrs['positive'] = 'down'

        dset = xr.Dataset(
            data_vars=data_vars,
            coords=coords,
        )

        return dset
