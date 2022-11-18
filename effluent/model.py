import io
import numpy as np
from scipy.integrate import solve_ivp
import xarray as xr
import logging
from pathlib import Path
import netCDF4 as nc


logger = logging.getLogger(__name__)


class Model:
    def __init__(self, fname_or_dict):
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


def load_config(fname_or_dict):
    """Load and parse input config, and convert to internal config format"""

    if isinstance(fname_or_dict, dict):
        input_conf = fname_or_dict
    else:
        import yaml
        with open(fname_or_dict, encoding='utf-8') as fp:
            input_conf = yaml.safe_load(fp)

    # noinspection PyDictCreation
    conf = {}

    conf['pipe'] = {}
    conf['pipe']['format'] = 'dict'
    conf['pipe']['time'] = input_conf['pipe']['time']
    conf['pipe']['flow'] = input_conf['pipe']['flow']
    conf['pipe']['dens'] = input_conf['pipe']['dens']
    conf['pipe']['decline'] = input_conf['pipe']['decline']
    conf['pipe']['diam'] = input_conf['pipe']['diam']
    conf['pipe']['depth'] = input_conf['pipe']['depth']

    conf['ambient'] = {}
    conf['ambient']['format'] = 'dict'
    conf['ambient']['time'] = input_conf['ambient']['time']
    conf['ambient']['depth'] = input_conf['ambient']['depth']
    conf['ambient']['coflow'] = input_conf['ambient']['coflow']
    conf['ambient']['crossflow'] = input_conf['ambient']['crossflow']
    conf['ambient']['dens'] = input_conf['ambient']['dens']

    conf['solver'] = {}
    conf['solver']['resolution'] = input_conf['output']['resolution']
    conf['solver']['stagnation'] = input_conf['output']['stagnation']

    conf['output'] = {}
    conf['output']['file'] = input_conf['output']['file']
    conf['output']['format'] = Path(input_conf['output']['file']).suffix[1:]

    conf['timestepper'] = {}
    conf['timestepper']['frequency'] = input_conf['output']['frequency']
    conf['timestepper']['stop'] = input_conf['output']['stop']

    return conf


class Pipe:
    def __init__(self, dset):
        self.dset = dset

    @staticmethod
    def from_config(conf):
        fmt = conf.pop('format')
        factories = {'dict': Pipe.from_mapping}
        factory = factories[fmt]
        return factory(**conf)

    @staticmethod
    def from_mapping(time, flow, dens, decline, diam, depth):
        theta = decline * np.pi / 180
        time = np.array(time)
        assert np.all(np.diff(time) > 0), "time values must be strictly increasing"
        shp = (len(time), )
        flow = np.broadcast_to(flow, shp)
        dens = np.broadcast_to(dens, shp)
        u = flow * np.cos(theta)
        w = flow * np.sin(theta)

        dset = xr.Dataset(
            data_vars=dict(
                u=xr.Variable('time', u),
                w=xr.Variable('time', w),
                dens=xr.Variable('time', dens),
                diam=xr.Variable((), diam),
            ),
            coords=dict(
                time=xr.Variable('time', time),
                depth=xr.Variable((), depth)
            ),
        )
        return Pipe(dset)

    def select(self, time):
        return interp_time(self.dset, time)


class Ambient:
    def __init__(self, dset):
        self.dset = dset

    @staticmethod
    def from_config(conf):
        fmt = conf.pop('format')
        factories = {'dict': Ambient.from_mapping}
        factory = factories[fmt]
        return factory(**conf)

    @staticmethod
    def from_mapping(time, depth, coflow, crossflow, dens):
        depth = np.array(depth)
        time = np.array(time)
        assert np.all(np.diff(time) > 0), "time values must be strictly increasing"
        shp = (len(time), len(depth))
        u = np.broadcast_to(coflow, shp)
        v = np.broadcast_to(crossflow, shp)

        dset = xr.Dataset(
            data_vars=dict(
                u=xr.Variable(('time', 'depth'), u),
                v=xr.Variable(('time', 'depth'), v),
                dens=xr.Variable(('time', 'depth'), dens),
            ),
            coords=dict(
                depth=xr.Variable('depth', depth),
                time=xr.Variable('time', time),
            )
        )
        return Ambient(dset)

    def select(self, time):
        return interp_time(self.dset, time)


def interp_time(dset, time):
    time_min = dset.time[0].values.item()
    time_max = dset.time[-1].values.item()
    clipped_time = np.clip(time, time_min, time_max)
    return dset.interp(time=clipped_time)


class Output:
    @staticmethod
    def from_config(conf):
        subclasses = {'csv': OutputCSV, 'nc': OutputNC}
        subclass = subclasses.get(conf.pop('format'), subclasses['nc'])
        return subclass(**conf)


class OutputCSV(Output):
    def __init__(self, file, diskless=False):
        self.file = file
        self.dset = None
        self._blank_file = True
        self.diskless = diskless

    def __enter__(self):
        if self.diskless:
            self.dset = io.StringIO()
        else:
            self.dset = open(self.file, 'w', encoding='utf-8', newline='\n')
        self._blank_file = True
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        if self.dset is not None:
            self.dset.close()
            self.dset = None

    def write(self, time, result):
        df = result.to_dataframe()
        df['release_time'] = time
        df = df.reset_index().set_index(['release_time', 't'])

        # Append result to file, write headers only if blank file
        df.to_csv(self.dset, line_terminator='\n', header=self._blank_file)
        self._blank_file = False


class OutputNC(Output):
    def __init__(self, file, diskless=False):
        self.file = file
        self.dset = None
        self._blank_file = True
        self.diskless = diskless

    def __enter__(self):
        self.dset = nc.Dataset(filename=self.file, mode='w', diskless=self.diskless)
        self._blank_file = True
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        if self.dset is not None:
            self.dset.close()
            self.dset = None

    def write(self, time, result):
        result = result.assign_coords(release_time=time)
        result = result.expand_dims('release_time')
        if self._blank_file:
            self._append_attributes(result)
            write_xr_to_nc(result, self.dset)
            self._blank_file = False
        else:
            append_xr_to_nc(result, self.dset)

    @staticmethod
    def _append_attributes(r):
        r.encoding['unlimited_dims'] = ['release_time']

        r['x'].attrs['long_name'] = 'centerline x coordinate'
        r['y'].attrs['long_name'] = 'centerline y coordinate'
        r['z'].attrs['long_name'] = 'centerline z coordinate'
        r['u'].attrs['long_name'] = 'velocity in x direction'
        r['v'].attrs['long_name'] = 'velocity in y direction'
        r['w'].attrs['long_name'] = 'velocity in z direction'
        r['density'].attrs['long_name'] = 'mass density of fluid'
        r['radius'].attrs['long_name'] = 'radius of top hat profile'
        r['t'].attrs['long_name'] = 'time since release'

        r['x'].attrs['units'] = 'm'
        r['y'].attrs['units'] = 'm'
        r['z'].attrs['units'] = 'm'
        r['u'].attrs['units'] = 'm/s'
        r['v'].attrs['units'] = 'm/s'
        r['w'].attrs['units'] = 'm/s'
        r['density'].attrs['units'] = 'kg/m^3'
        r['radius'].attrs['units'] = 'm'
        r['t'].attrs['units'] = 's'

        r['u'].attrs['standard_name'] = 'sea_water_x_velocity'
        r['v'].attrs['standard_name'] = 'sea_water_x_velocity'
        r['w'].attrs['standard_name'] = 'downward_sea_water_velocity'
        r['z'].attrs['standard_name'] = 'depth'
        r['density'].attrs['standard_name'] = 'sea_water_density'

        r['w'].attrs['positive'] = 'down'
        r['z'].attrs['positive'] = 'down'

        r.coords['release_time'].attrs['long_name'] = 'time of release'
        r.coords['release_time'].attrs['units'] = 's'


def write_xr_to_nc(xr_dset: xr.Dataset, nc_dset: nc.Dataset):
    unlimited_dims = xr_dset.encoding.get('unlimited_dims', [])

    # Write dimensions
    for name, size in xr_dset.dims.items():
        if name in unlimited_dims:
            size = None
        nc_dset.createDimension(name, size)

    # Write variables
    for name, xr_var in xr_dset.variables.items():
        nc_var = nc_dset.createVariable(
            varname=name,
            datatype=xr_var.dtype,
            dimensions=xr_var.dims,
            fill_value=False,
        )
        nc_var[:] = xr_var.values
        nc_var.setncatts(xr_var.attrs)

    # Write dataset attributes
    nc_dset.setncatts(xr_dset.attrs)


def append_xr_to_nc(xr_dset: xr.Dataset, nc_dset: nc.Dataset):
    unlim_dims = [k for k, v in nc_dset.dimensions.items() if v.isunlimited()]
    unlim_dim = unlim_dims[0] if len(unlim_dims) > 0 else None
    unlim_vars = [k for k, v in xr_dset.variables.items() if v.dims[0] == unlim_dim]

    num_old_items = nc_dset.dimensions[unlim_dim].size
    num_new_items = xr_dset.dims.get(unlim_dim, 0)
    num_items = num_old_items + num_new_items

    # Append data
    for name in unlim_vars:
        xr_var = xr_dset[name]
        nc_var = nc_dset.variables[name]
        nc_var[num_old_items:num_items] = xr_var.values


class Solver:
    def __init__(self, resolution, stagnation):
        self.resolution = resolution
        self.stagnation = stagnation

    @staticmethod
    def from_config(conf):
        return Solver(**conf)

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
        res_t, res_y = result.t, result.y

        # Organize result
        dset = xr.Dataset(
            data_vars={v: xr.Variable('t', res_y[i]) for i, v in enumerate(varnames)},
            coords=dict(t=res_t),
        )

        return dset
