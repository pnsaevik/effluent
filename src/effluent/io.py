import abc

import netCDF4 as nc
import numpy as np
import pandas as pd
import tomli as toml
import xarray as xr


def load_config(fname_or_dict):
    """Load and parse input config, and convert to internal config format"""

    if isinstance(fname_or_dict, dict):
        input_conf = fname_or_dict
    else:
        with open(fname_or_dict, 'rb') as fp:
            input_conf = toml.load(fp)

    # noinspection PyDictCreation
    conf = {}
    conf['pipe'] = input_conf['pipe']
    conf['ambient'] = input_conf['ambient']
    conf['output'] = input_conf['output']

    # --- Solver ---

    conf['solver'] = {}
    # Copy a selection of output parameters
    for k in ['stagnation', 'resolution']:
        conf['solver'][k] = input_conf['output'][k]
    # Copy model parameters
    for k, v in input_conf.get('model', {}).items():
        conf['solver'][k] = v
    # Copy solver parameters
    for k, v in input_conf.get('solver', {}).items():
        conf['solver'][k] = v

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
        elif 'nc' in conf:
            return Pipe.from_nc_file(**conf['nc'])
        else:
            return Pipe.from_mapping(**conf)

    @staticmethod
    def from_nc_file(file):
        dset = xr.load_dataset(file)
        return Pipe.from_dataset(dset)

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
        return Pipe.from_dataframe(df)

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
        self._dset = dset
        self._tmin = dset.time[0].values.item()
        self._tmax = dset.time[-1].values.item()

    @staticmethod
    def from_config(conf):
        if 'csv' in conf:
            return Ambient.from_csv_file(**conf['csv'])
        elif 'nc' in conf:
            return Ambient.from_nc_file(**conf['nc'])
        else:
            return Ambient.from_mapping(**conf)

    @staticmethod
    def from_nc_file(file):
        dset = xr.load_dataset(file)
        return Ambient.from_dataset(dset)

    @staticmethod
    def from_dataframe(df):
        dset = xr.Dataset.from_dataframe(df)
        return Ambient.from_dataset(dset)

    @staticmethod
    def from_dataset(dset):
        dset = dset.rename_vars(coflow='u', crossflow='v')
        time = dset.time.values
        assert np.all(np.diff(time) > 0), "time values must be strictly increasing"
        return Ambient(dset)

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
        return Ambient.from_dataframe(df)

    @staticmethod
    def from_mapping(time, depth, coflow, crossflow, dens):
        shp = (len(time), len(depth))
        u = np.broadcast_to(coflow, shp)
        v = np.broadcast_to(crossflow, shp)
        d = np.broadcast_to(dens, shp)

        dset = xr.Dataset(
            coords=dict(time=time, depth=depth),
            data_vars=dict(
                coflow=xr.Variable(('time', 'depth'), u),
                crossflow=xr.Variable(('time', 'depth'), v),
                dens=xr.Variable(('time', 'depth'), d),
            ),
        )
        return Ambient.from_dataset(dset)

    def select(self, time):
        clipped_time = np.clip(time, self._tmin, self._tmax)
        return self._dset.interp(time=clipped_time)


class Output:
    @staticmethod
    def from_config(conf):
        if 'csv' in conf:
            return OutputCSV.from_config(conf)
        elif 'nc' in conf:
            return OutputNC.from_config(conf)
        else:
            raise ValueError("No output file name given")

    @abc.abstractmethod
    def write(self, time, result):
        return NotImplementedError


class OutputCSV(Output):
    def __init__(self, file, variables):
        self.variables = variables
        self._blank_file = True

        if isinstance(file, str):
            self.file = file
            self.dset = None

        elif hasattr(file, 'write') and callable(file.write):
            # Diskless mode
            self.file = None
            self.dset = file

        else:
            raise TypeError(f'Expected file name or stream, found "{type(file)}"')

    def __enter__(self):
        if self.dset is None:
            self.dset = open(self.file, 'w', encoding='utf-8', newline='\n')
            self._blank_file = True
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        if self.dset is not None:
            self.dset.close()
            self.dset = None

    @staticmethod
    def from_config(conf):
        out = OutputCSV(
            file=conf['csv']['file'],
            variables=conf.get('variables', None)
        )
        return out

    def write(self, time, result):
        df = result.to_dataframe()
        df['release_time'] = time
        df = df.reset_index().set_index(['release_time', 't']).reset_index()

        # Drop non-output variables
        if self.variables is not None:
            df = df[self.variables]

        # Append result to file, write headers only if blank file
        df.to_csv(
            self.dset,
            line_terminator='\n',
            header=self._blank_file,
            float_format='%.10g',
            index=False,
        )
        self._blank_file = False


class OutputNC(Output):
    def __init__(self, file, variables):
        self.variables = variables
        self.dset = None
        self._blank_file = True

        if isinstance(file, str):
            self.fname = file
            self.diskless = False
            self.xr_dset = None

        elif isinstance(file, xr.Dataset):
            # Diskless mode: Data is written to memory, and to an xarray.Dataset on exit
            from uuid import uuid4
            self.fname = uuid4()
            self.diskless = True
            self.xr_dset = file

        else:
            raise TypeError(f'Unknown file type: {type(file)}')

    def __enter__(self):
        self.dset = nc.Dataset(filename=self.fname, mode='w', diskless=self.diskless)
        self._blank_file = True
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        if self.dset is not None:
            if self.xr_dset is not None:
                write_nc_to_xr(self.dset, self.xr_dset)

            self.dset.close()
            self.dset = None

    @staticmethod
    def from_config(conf):
        out = OutputNC(
            file=conf['nc']['file'],
            variables=conf.get('variables', None)
        )
        return out

    def write(self, time, result):
        result = result.assign_coords(release_time=time)
        result = result.expand_dims('release_time')

        if self._blank_file:
            self._append_attributes(result)

        # Drop non-output variables
        if self.variables is not None:
            non_output_vars = [v for v in result.variables if v not in self.variables]
            result = result.drop_vars(non_output_vars)

        if self._blank_file:
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


def write_nc_to_xr(nc_dset: nc.Dataset, xr_dset: xr.Dataset):
    # Write variables
    for name, nc_var in nc_dset.variables.items():
        xr_var = xr.Variable(
            dims=nc_var.dimensions,
            data=nc_var[:],
            attrs={k: nc_var.getncattr(k) for k in nc_var.ncattrs()},
        )
        xr_dset[name] = xr_var

    # Write dataset attributes
    for k in nc_dset.ncattrs():
        xr_dset.attrs[k] = nc_dset.getncattr(k)


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
