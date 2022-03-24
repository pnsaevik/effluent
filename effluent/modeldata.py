import xarray as xr
import numpy as np


def _create_profile(z, dens, u=None, v=None, w=None, lat=None, lon=None, date=None):
    u = u or np.zeros(len(z))
    v = v or np.zeros(len(z))
    w = w or np.zeros(len(z))
    date = date or np.array(np.nan)
    lat = lat or np.nan
    lon = lon or np.nan

    return xr.Dataset(
        data_vars=dict(
            u=xr.Variable(
                dims='z',
                data=u,
                attrs=dict(
                    standard_name='eastward_sea_water_velocity',
                    units='m/s',
                ),
            ),
            v=xr.Variable(
                dims='z',
                data=v,
                attrs=dict(
                    standard_name='northward_sea_water_velocity',
                    units='m/s',
                ),
            ),
            w=xr.Variable(
                dims='z',
                data=w,
                attrs=dict(
                    standard_name='upward_sea_water_velocity',
                    units='m/s',
                ),
            ),
            dens=xr.Variable(
                dims='z',
                data=dens,
                attrs=dict(
                    standard_name='sea_water_density',
                    units='kg m-3',
                ),
            ),
        ),
        coords=dict(
            z=xr.Variable(
                dims='z',
                data=z,
                attrs=dict(
                    standard_name='depth',
                    units='m',
                )
            ),
            lon=xr.Variable(
                dims=(),
                data=lon,
                attrs=dict(
                    standard_name='longitude',
                    units='degrees_east',
                ),
            ),
            lat=xr.Variable(
                dims=(),
                data=lat,
                attrs=dict(
                    standard_name='latitude',
                    units='degrees_north',
                ),
            ),
            date=xr.Variable(
                dims=(),
                data=date.astype('datetime64[s]'),
                attrs=dict(
                    standard_name='time',
                    calendar='standard',
                    units='seconds since 1970-01-01',
                ),
            ),
        ),
    )


class Modeldata:
    def __init__(self):
        pass

    @staticmethod
    def from_profile(z, dens, u=None, v=None, w=None, lat=None, lon=None, date=None):
        return ModeldataSingleProfile(z, dens, u, v, w, lat, lon, date)

    @staticmethod
    def from_roms(file_names):
        return ModeldataRoms(file_names)

    def profile(self, lat, lon, date):
        raise NotImplementedError("Use from_* methods to generate an instance")


# --- ROMS functions ---

class ModeldataRoms(Modeldata):
    def __init__(self, file_names):
        super().__init__()
        self.filetab = _filetab_create(file_names)
        self.horizontal_transform = _roms_horizontal_transform(file_names[0])

    def profile(self, lat, lon, date):

        pass

    def _from_date_to_index(self, date):
        dates = self.filetab.time_val.values
        i = np.searchsorted(dates, date)
        row = self.filetab.iloc[i]
        return row.file_index, row.time_index


def _roms_horizontal_transform(file_name):
    from .utils import bilin_inv

    with xr.open_dataset(file_name) as dset:
        lat_rho = dset.variables['lat_rho'].values
        lon_rho = dset.variables['lon_rho'].values

    previous_latlon = (None, None)
    previous_xy = [None, None]

    def latlon_to_xy(lat, lon):
        if previous_latlon != (lat, lon):
            x, y = bilin_inv(lat, lon, lat_rho, lon_rho)
            previous_xy[0] = x
            previous_xy[1] = y
        return previous_xy

    return latlon_to_xy


def _filetab_create(file_names):
    if isinstance(file_names, str):
        file_names = [file_names]

    import pandas as pd

    file_index_series = []
    file_name_series = []
    time_index_series = []
    time_val_series = []
    for i, fname in enumerate(file_names):
        with xr.open_dataset(fname) as dset:
            t = dset.variables['ocean_time'].values.astype('datetime64[us]').tolist()
        time_val_series += t
        time_index_series += list(range(len(t)))
        file_name_series += [fname] * len(t)
        file_index_series += [i] * len(t)

    df = pd.DataFrame(
        dict(
            file_name=file_name_series,
            file_index=file_index_series,
            time_val=time_val_series,
            time_index=time_index_series,
        )
    )

    df.sort_values(by='time_val', inplace=True)
    df.drop_duplicates(subset='time_val', inplace=True)

    return df


# --- End ROMS functions


class ModeldataSingleProfile(Modeldata):
    def __init__(self, z, dens, u=None, v=None, w=None, lat=None, lon=None, date=None):
        super().__init__()
        self._profile = _create_profile(z, dens, u, v, w, lat, lon, date)

    def profile(self, lat=None, lon=None, date=None):
        return self._profile
