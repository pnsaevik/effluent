import io
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import xarray as xr
import logging
from pathlib import Path
import netCDF4 as nc
import tomli as toml


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
        with open(fname_or_dict, 'rb') as fp:
            input_conf = toml.load(fp)

    # noinspection PyDictCreation
    conf = {}

    # --- Pipe ---

    c = input_conf['pipe']

    # Interpret input as a file name if there is only one input string
    if isinstance(c, str):
        c = dict(
            file=c,
            format=Path(c).suffix[1:]
        )

    # Add standard format specifier if not existing
    if 'format' not in c:
        c['format'] = 'dict'

    conf['pipe'] = c

    # --- Ambient ---

    c = input_conf['ambient']

    # Interpret input as a file name if there is only one input string
    if isinstance(c, str):
        c = dict(
            file=c,
            format=Path(c).suffix[1:]
        )

    # Add standard format specifier if not existing
    if 'format' not in c:
        c['format'] = 'dict'

    conf['ambient'] = c

    # --- Solver ---

    stagnation = input_conf['output']['stagnation']
    resolution = input_conf['output']['resolution']

    conf['solver'] = {}
    conf['solver']['steps'] = np.arange(0, stagnation, resolution)

    # --- Output ---

    conf['output'] = {}
    conf['output']['file'] = input_conf['output']['file']
    conf['output']['format'] = Path(input_conf['output']['file']).suffix[1:]

    # --- Time stepper ---

    conf['timestepper'] = {}
    conf['timestepper']['frequency'] = input_conf['output']['frequency']
    conf['timestepper']['stop'] = input_conf['output']['stop']

    return conf


class Pipe:
    def __init__(self, dset):
        self._dset = dset
        self._time_min = dset.time[0].values.item()
        self._time_max = dset.time[-1].values.item()

    @staticmethod
    def from_config(conf):
        if 'csv' in conf:
            return Pipe.from_csv_file(**conf['csv'])
        else:
            return Pipe.from_mapping(**conf)

    @staticmethod
    def from_csv_file(file):
        df = pd.read_csv(
            file,
            sep=',',
            index_col='time',
            header=0,
            skipinitialspace=True,
            skip_blank_lines=True,
            comment='#',
        )

        df['u'], df['w'] = Pipe._compute_uw(df['flow'].values, df['decline'].values)
        df = df.drop(columns=['flow', 'decline'])
        dset = xr.Dataset.from_dataframe(df)
        return Pipe(dset)

    @staticmethod
    def _compute_uw(flow, decline):
        theta = decline * np.pi / 180
        u = flow * np.cos(theta)
        w = flow * np.sin(theta)
        return u, w

    @staticmethod
    def from_dataframe(df):
        dset = xr.Dataset.from_dataframe(df)
        return Pipe.from_dataset(dset)

    @staticmethod
    def from_dataset(dset):
        time = dset.time.values
        assert np.all(np.diff(time) > 0), "time values must be strictly increasing"
        u, w = Pipe._compute_uw(dset.flow.values, dset.decline.values)
        dset['u'] = xr.Variable('time', u)
        dset['w'] = xr.Variable('time', w)
        return Pipe(dset)

    @staticmethod
    def from_mapping(time, flow, dens, decline, diam, depth):
        index = pd.Index(data=time, name='time')
        data = dict(flow=flow, dens=dens, diam=diam, depth=depth, decline=decline)
        df = pd.DataFrame(data, index=index)
        return Pipe.from_dataframe(df)

    def select(self, time):
        clipped_time = np.clip(time, self._time_min, self._time_max)
        return self._dset.interp(time=clipped_time)


class Ambient:
    def __init__(self, dset):
        self.dset = dset

    @staticmethod
    def from_config(conf):
        fmt = conf.pop('format')
        factories = {'dict': Ambient.from_mapping, 'csv': Ambient.from_csv_file}
        factory = factories[fmt]
        return factory(**conf)

    @staticmethod
    def from_csv_file(file):
        df = pd.read_csv(
            file,
            sep=',',
            index_col=('time', 'depth'),
            header=0,
            skipinitialspace=True,
            skip_blank_lines=True,
            comment='#',
        )

        dset = xr.Dataset.from_dataframe(df)
        dset = dset.rename_vars(coflow='u', crossflow='v')
        return Ambient(dset)

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
        df.to_csv(
            self.dset,
            line_terminator='\n',
            header=self._blank_file,
            float_format='%.10g',
        )
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
        z_min, z_max = self.ambient.depth[[0, -1]].values
        ambient = self.ambient.interp(depth=np.clip(z, z_min, z_max))
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
