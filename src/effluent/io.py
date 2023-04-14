import abc

import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr


class Pipe:
    def __init__(self, dset):
        self._dset = dset
        self._time_min = dset.time[0].values
        self._time_max = dset.time[-1].values

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
        df = read_csv(file)
        return Pipe.from_dataframe(df)

    @staticmethod
    def _compute_uw(flow, decline, diam):
        area = np.pi * diam * diam * 0.25
        speed = flow / area
        theta = decline * np.pi / 180
        u = speed * np.cos(theta)
        w = speed * np.sin(theta)
        return u, w

    @staticmethod
    def from_dataframe(df):
        df['time'] = df['time'].values.astype('datetime64')
        df = df.set_index('time')
        dset = xr.Dataset.from_dataframe(df)
        return Pipe.from_dataset(dset)

    @staticmethod
    def from_dataset(dset):
        u, w = Pipe._compute_uw(dset.flow.values, dset.decline.values, dset.diam.values)
        dset['u'] = xr.Variable('time', u)
        dset['w'] = xr.Variable('time', w)

        # Compute density from temp and salt if not present
        if 'dens' not in dset:
            from . import eos
            dens = eos.roms_rho(
                temp=dset['temp'].values,
                salt=dset['salt'].values,
                depth=dset['depth'].values,
            )
            dset = dset.assign(dens=xr.Variable(dset['temp'].dims, dens))
        return Pipe(dset)

    @staticmethod
    def from_mapping(time, flow, decline, diam, depth, dens=None, salt=None, temp=None):
        data = dict(time=time, flow=flow, diam=diam, depth=depth, decline=decline)
        if salt is not None:
            data['salt'] = salt
        if temp is not None:
            data['temp'] = temp
        if dens is not None:
            data['dens'] = dens
        df = pd.DataFrame(data)
        return Pipe.from_dataframe(df)

    def select(self, time):
        if self._dset.dims['time'] == 1:
            # No interpolation is possible if there is only 1 time entry
            return self._dset.isel(time=0)

        clipped_time = np.clip(time, self._time_min, self._time_max)
        return self._dset.interp(time=clipped_time)


def read_csv(file):
    return pd.read_csv(
        file,
        sep=',',
        header=0,
        skipinitialspace=True,
        skip_blank_lines=True,
        comment='#',
        converters=dict(time=np.datetime64),
    )


class Ambient:
    @abc.abstractmethod
    def select(self, time):
        return NotImplementedError

    @staticmethod
    def from_config(conf):
        if 'csv' in conf:
            return Ambient.from_csv_file(**conf['csv'])
        elif 'nc' in conf:
            return Ambient.from_nc_file(**conf['nc'])
        elif 'roms' in conf:
            return AmbientRoms(**conf['roms'])
        else:
            return Ambient.from_mapping(**conf)

    @staticmethod
    def from_nc_file(file):
        dset = xr.load_dataset(file)
        return Ambient.from_dataset(dset)

    @staticmethod
    def from_dataframe(df):
        df['time'] = df['time'].values.astype('datetime64')
        df = df.set_index(['time', 'depth'])

        dset = xr.Dataset.from_dataframe(df)
        return Ambient.from_dataset(dset)

    @staticmethod
    def from_dataset(dset):
        dset = dset.rename_vars(coflow='u', crossflow='v')
        time = dset.time.values
        assert np.all(np.diff(time).astype('int64') > 0), "time values must be strictly increasing"

        # Compute density from temp and salt if not present
        if 'dens' not in dset:
            from . import eos
            data = eos.roms_rho(
                temp=dset['temp'].values,
                salt=dset['salt'].values,
                depth=dset['depth'].values,
            )
            dset = dset.assign(dens=xr.Variable(dset['temp'].dims, data))

        return AmbientXarray(dset)

    @staticmethod
    def from_csv_file(file):
        df = read_csv(file)
        return Ambient.from_dataframe(df)

    @staticmethod
    def from_mapping(time, depth, coflow, crossflow, dens=None, salt=None, temp=None):
        shp = (len(time), len(depth))

        variables = dict(coflow=coflow, crossflow=crossflow, dens=dens, salt=salt,
                         temp=temp)

        dset = xr.Dataset(coords=dict(time=time, depth=depth))
        for k, v in variables.items():
            if v is None:
                continue
            data = np.broadcast_to(v, shp)
            dset = dset.assign(**{k: xr.Variable(('time', 'depth'), data)})

        return Ambient.from_dataset(dset)

    def close(self):
        pass


class AmbientXarray(Ambient):
    def __init__(self, dset):
        self._dset = dset
        self._tmin = dset.time[0].values
        self._tmax = dset.time[-1].values

    def select(self, time):
        if self._dset.dims['time'] == 1:
            # No interpolation is possible if there is only 1 time entry
            return self._dset.isel(time=0)

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

    def close(self):
        pass


class OutputCSV(Output):
    def __init__(self, file, variables, float_format, separator):
        self.variables = variables
        self.float_format = float_format
        self.separator = separator
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
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def open(self):
        if self.dset is None:
            self.dset = open(self.file, 'w', encoding='utf-8', newline='\n')
            self._blank_file = True

    def close(self):
        if self.dset is not None:
            self.dset.close()
            self.dset = None

    @staticmethod
    def from_config(conf):
        out = OutputCSV(
            file=conf['csv']['file'],
            variables=conf.get('variables', None),
            float_format=conf['csv'].get('float_format', '%.10g'),
            separator=conf.get('separator', ','),
        )
        return out

    def write(self, time, result):
        self.open()  # Lazy opening: Only effective if first time

        df = result.to_dataframe()
        df['release_time'] = time
        df = df.reset_index().set_index(['release_time', 't']).reset_index()

        # Drop non-output variables
        if self.variables is not None:
            df = df[self.variables]

        # Append result to file, write headers only if blank file
        df.to_csv(
            self.dset,
            lineterminator='\n',
            header=self._blank_file,
            float_format=self.float_format,
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
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def open(self):
        if self.dset is None:
            self.dset = nc.Dataset(filename=self.fname, mode='w', diskless=self.diskless)
            self._blank_file = True

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
        self.open()  # Lazy opening: Only effective if first time

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


class AmbientRoms(Ambient):
    def __init__(self, file, latitude, longitude, azimuth):
        self.file = file
        self.latitude = latitude
        self.longitude = longitude
        self.azimuth = azimuth

        self.dset = None
        self._tmin = None
        self._tmax = None

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def open(self):
        if self.dset is None:
            import effluent.roms
            dset = effluent.roms.open_location(
                file=self.file,
                latitude=self.latitude,
                longitude=self.longitude,
                azimuth=self.azimuth,
            )

            keep_vars = ['time', 'depth', 'u', 'v', 'dens']
            drop_vars = [v for v in dset.variables if v not in keep_vars]
            dset = dset.drop_vars(drop_vars)

            self.dset = dset
            self._tmin = dset.time[0].values
            self._tmax = dset.time[-1].values

    def close(self):
        if self.dset is not None:
            self.dset.close()
            self.dset = None

    def select(self, time):
        self.open()

        clipped_time = np.clip(time, self._tmin, self._tmax)
        return self.dset.interp(time=clipped_time)
