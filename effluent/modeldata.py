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
    def from_single_profile(z, dens, u=None, v=None, w=None, lat=None, lon=None, date=None):
        return ModeldataSingleProfile(z, dens, u, v, w, lat, lon, date)

    def profile(self, lat, lon, date):
        raise NotImplementedError("Use from_* methods to generate an instance")


class ModeldataSingleProfile(Modeldata):
    def __init__(self, z, dens, u=None, v=None, w=None, lat=None, lon=None, date=None):
        super().__init__()
        self._profile = _create_profile(z, dens, u, v, w, lat, lon, date)

    def profile(self, lat=None, lon=None, date=None):
        return self._profile
