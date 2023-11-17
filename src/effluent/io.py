"""
The package contains functions and classes for reading and writing simulation data.
"""

import abc
import numpy as np
import xarray as xr
import pandas as pd
import effluent.eos
import netCDF4 as nc
import uuid
import effluent.roms
import cftime


class Pipe:
    """
    Data about the pipe and the effluent release.

    The constructor takes an xarray.Dataset object as input. If the source data is
    in the form of a data file, the factory method should be used instead.

    :param dset: An xarray.Dataset object with variables ``depth``, ``u``,
        ``w``, ``dens`` and ``diam``, all indexed by the coordinate ``time``.
    """

    def __init__(self, dset):
        self._dset = dset
        self._time_min = dset.time[0].values
        self._time_max = dset.time[-1].values

    @staticmethod
    def from_config(conf) -> "Pipe":
        """
        Initialize using :doc:`configuration parameters </config/pipe>`

        :param conf: A dict of configuration parameters
        :return: An initialized object
        """
        if 'csv' in conf:
            return Pipe.from_csv_file(**conf['csv'])
        elif 'nc' in conf:
            return Pipe.from_nc_file(**conf['nc'])
        else:
            return Pipe.from_mapping(**conf)

    @staticmethod
    def from_nc_file(file) -> "Pipe":
        dset = xr.load_dataset(file)
        return Pipe.from_dataset(dset)

    @staticmethod
    def from_csv_file(file) -> "Pipe":
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
    def from_dataframe(df) -> "Pipe":
        df['time'] = df['time'].values.astype('datetime64')
        df = df.set_index('time')
        dset = xr.Dataset.from_dataframe(df)
        return Pipe.from_dataset(dset)

    @staticmethod
    def from_dataset(dset) -> "Pipe":
        u, w = Pipe._compute_uw(dset.flow.values, dset.decline.values, dset.diam.values)
        dset['u'] = xr.Variable('time', u)
        dset['w'] = xr.Variable('time', w)

        # Compute density from temp and salt if not present
        if 'dens' not in dset:
            dens = effluent.eos.roms_rho(
                temp=dset['temp'].values,
                salt=dset['salt'].values,
                depth=dset['depth'].values,
            )
            dset = dset.assign(dens=xr.Variable(dset['temp'].dims, dens))
        return Pipe(dset)

    @staticmethod
    def from_mapping(time, flow, decline, diam, depth, dens=None, salt=None, temp=None) -> "Pipe":
        data = dict(time=time, flow=flow, diam=diam, depth=depth, decline=decline)
        if salt is not None:
            data['salt'] = salt
        if temp is not None:
            data['temp'] = temp
        if dens is not None:
            data['dens'] = dens
        df = pd.DataFrame(data)
        return Pipe.from_dataframe(df)

    def select(self, time) -> xr.Dataset:
        """
        Interpolate pipe parameters to a specific point in time

        :param time: A time in numpy.datetime64 format
        :return: An xarray.Dataset object with variables ``depth``, ``u``, ``w``,
            ``dens`` and ``diam``.
        """

        if self._dset.dims['time'] == 1:
            # No interpolation is possible if there is only 1 time entry
            return self._dset.isel(time=0)

        clipped_time = np.clip(time, self._time_min, self._time_max)
        return self._dset.interp(time=clipped_time)


def read_csv(file) -> pd.DataFrame:
    """
    Read csv file, and return a pandas.DataFrame

    :param file: File name
    :return: A pandas.DataFrame object
    """
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
    """
    Data about the ambient ocean.

    This is an abstract base class with no explicit constructor. To initialize an
    instance of the class, use the factory method.
    """

    @abc.abstractmethod
    def select(self, time) -> xr.Dataset:
        """
        Compute the ambient conditions at a specific time.

        :param time: A numpy.datetime64 object
        :return: An xarray.Dataset object with variables ``u``, ``v`` and ``dens``, all
            indexed by the coordinate ``depth``.
        """
        # noinspection PyTypeChecker
        return NotImplementedError

    @staticmethod
    def from_config(conf) -> "Ambient":
        """
        Initialize using :doc:`configuration parameters </config/ambient>`

        :param conf: A dict of configuration parameters
        :return: An initialized object
        """

        if 'csv' in conf:
            return Ambient.from_csv_file(**conf['csv'])
        elif 'nc' in conf:
            return Ambient.from_nc_file(**conf['nc'])
        elif 'roms' in conf:
            return AmbientRoms(**conf['roms'])
        else:
            return Ambient.from_mapping(**conf)

    @staticmethod
    def from_nc_file(file) -> "Ambient":
        dset = xr.load_dataset(file)
        return Ambient.from_dataset(dset)

    @staticmethod
    def from_dataframe(df) -> "Ambient":
        df['time'] = df['time'].values.astype('datetime64')
        df = df.set_index(['time', 'depth'])

        dset = xr.Dataset.from_dataframe(df)
        return Ambient.from_dataset(dset)

    @staticmethod
    def from_dataset(dset) -> "Ambient":
        dset = dset.rename_vars(coflow='u', crossflow='v')
        time = dset.time.values
        assert np.all(np.diff(time).astype('int64') > 0), "time values must be strictly increasing"

        # Compute density from temp and salt if not present
        if 'dens' not in dset:
            data = effluent.eos.roms_rho(
                temp=dset['temp'].values,
                salt=dset['salt'].values,
                depth=dset['depth'].values,
            )
            dset = dset.assign(dens=xr.Variable(dset['temp'].dims, data))

        return AmbientXarray(dset)

    @staticmethod
    def from_csv_file(file) -> "Ambient":
        df = read_csv(file)
        return Ambient.from_dataframe(df)

    @staticmethod
    def from_mapping(time, depth, coflow, crossflow, dens=None, salt=None, temp=None) -> "Ambient":
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
        """
        Close the underlying data source
        """
        pass


class AmbientXarray(Ambient):
    """
    Data about the ambient ocean, from in-memory dataset.

    This subclass is using an xarray.Dataset object as its data source. The dataset
    should have variables ``u``, ``v`` and ``dens``, all indexed by the coordinates
    ``depth`` and ``time``.

    :param dset: An xarray.Dataset object
    """

    def __init__(self, dset):
        self._dset = dset
        self._tmin = dset.time[0].values
        self._tmax = dset.time[-1].values

    def select(self, time) -> xr.Dataset:
        if self._dset.dims['time'] == 1:
            # No interpolation is possible if there is only 1 time entry
            return self._dset.isel(time=0)

        clipped_time = np.clip(time, self._tmin, self._tmax)
        return self._dset.interp(time=clipped_time)


class Output:
    """
    Class for writing simulation output to disk

    This is an abstract base class with no explicit constructor. To initialize an
    instance of the class, use the factory method.
    """

    @staticmethod
    def from_config(conf) -> "Output":
        """
        Initialize using :doc:`configuration parameters </config/output>`

        :param conf: A dict of configuration parameters
        :return: An initialized object
        """

        if 'csv' in conf:
            return OutputCSV.from_config(conf)
        elif 'nc' in conf:
            return OutputNC.from_config(conf)
        else:
            raise ValueError("No output file name given")

    @abc.abstractmethod
    def write(self, time, result):
        """
        Write simulation results to disk

        :param time: Discharge time, as numpy.datetime64 object
        :param result: Simulation result, as returned by :func:`effluent.solver.Solver.solve`
        """
        return NotImplementedError

    def close(self):
        """
        Close the underlying data stream
        """
        pass


class OutputCSV(Output):
    """
    Class for writing simulation output to CSV file.

    The output file is created lazily upon the first write statement.

    :param file: Name of output file
    :param variables: A list of variable names to include
    :param float_format: Output format for float numbers
    :param separator: Symbol used as data separator
    """
    def __init__(self, file, variables=None, float_format="%.10g", separator=","):
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
    def from_config(conf) -> "OutputCSV":
        params = dict(file=conf['csv']['file'])
        if 'variables' in conf:
            params['variables'] = conf['variables']
        if 'float_format' in conf['csv']:
            params['float_format'] = conf['csv']['float_format']
        if 'separator' in conf['csv']:
            params['separator'] = conf['csv']['separator']

        return OutputCSV(**params)

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
            header=self._blank_file,
            float_format=self.float_format,
            index=False,
        )
        self._blank_file = False


class OutputNC(Output):
    """
    Class for writing simulation output to netCDF file.

    The output file is created lazily upon the first write statement.

    :param file: Name of output file
    :param variables: A list of variable names to include
    """

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
            self.fname = uuid.uuid4()
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
    def from_config(conf) -> "OutputNC":
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


def convert_to_nc_date(
        xr_var: xr.Variable, units: str = 'seconds since 1970-01-01',
        calendar: str = 'proleptic_gregorian', dtype='i8',
) -> xr.Variable:
    """
    Convert an xarray date variable to a CF-style date variable

    On default, the xarray date will be converted to seconds (int64) since the unix epoch
    (1970-01-01) using the proleptic gregorian calendar.

    :param xr_var: Input variable, with values of type numpy.datetime64
    :param units: Output units, as defined by CF conventions
    :param calendar: Output calendar, as defined by CF conventions
    :param dtype: Output data type, as defined by numpy conventions
    :return: Output variable, with values converted to seconds since epoch
    """
    dates = xr_var.values.astype('datetime64[us]').tolist()
    cf_dates = cftime.date2num(dates=dates, units=units, calendar=calendar).astype(dtype)
    new_attrs = {**xr_var.attrs, **dict(units=units, calendar=calendar)}
    return xr.Variable(dims=xr_var.dims, data=cf_dates, attrs=new_attrs)


def write_xr_to_nc(xr_dset: xr.Dataset, nc_dset: nc.Dataset):
    """
    Write data from an xarray.Dataset to a netCDF4.Dataset

    :param xr_dset: Input dataset
    :param nc_dset: Output dataset
    """
    unlimited_dims = xr_dset.encoding.get('unlimited_dims', [])

    # Write dimensions
    for name, size in xr_dset.dims.items():
        if name in unlimited_dims:
            size = None
        nc_dset.createDimension(name, size)

    # Write variables
    for name, xr_var in xr_dset.variables.items():
        if np.issubdtype(xr_var.dtype, np.datetime64):
            xr_var = convert_to_nc_date(xr_var)

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
    """
    Write data from a netCDF4.Dataset to a an xarray.Dataset

    :param nc_dset: Input dataset
    :param xr_dset: Output dataset
    """

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
    """
    Append data from an xarray.Dataset to a netCDF4.Dataset

    This method does not create new variables in the destination dataset, but only
    appends to the existing variables.

    :param xr_dset: Input dataset
    :param nc_dset: Output dataset
    """
    unlim_dims = [k for k, v in nc_dset.dimensions.items() if v.isunlimited()]
    unlim_dim = unlim_dims[0] if len(unlim_dims) > 0 else None
    unlim_vars = [k for k, v in xr_dset.variables.items() if v.dims[0] == unlim_dim]

    num_old_items = nc_dset.dimensions[unlim_dim].size
    num_new_items = xr_dset.dims.get(unlim_dim, 0)
    num_items = num_old_items + num_new_items

    # Append data
    for name in unlim_vars:
        xr_var = xr_dset.variables[name]
        if np.issubdtype(xr_var.dtype, np.datetime64):
            xr_var = convert_to_nc_date(xr_var, units=nc_dset[name].units,
                                        calendar=nc_dset[name].calendar)
        nc_var = nc_dset.variables[name]
        nc_var[num_old_items:num_items] = xr_var.values


class AmbientRoms(Ambient):
    """
    Data about the ambient ocean, from ROMS.

    This subclass is using ROMS output files as its data source. The constructor also
    needs the pipe position and azimuth, as described in the
    :doc:`documentation </config/ambient>` of the configuration file.

    The class lazily opens the underlying data source the first time it's needed.

    :param file: A set of ROMS files, specified with a wildcard string
    :param latitude: The pipe latitude
    :param longitude: The pipe longitude
    :param azimuth: The pipe azimuth (north = 0, east = 90)
    """

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
            dset = effluent.roms.load_location(
                file=self.file,
                lat=self.latitude,
                lon=self.longitude,
                az=self.azimuth,
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

    def select(self, time) -> xr.Dataset:
        self.open()

        clipped_time = np.clip(time, self._tmin, self._tmax)
        return self.dset.interp(time=clipped_time)
