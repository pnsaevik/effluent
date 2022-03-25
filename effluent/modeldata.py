import contextlib

import xarray as xr
import numpy as np


def _create_profile(z, dens, u=None, v=None, w=None, lat=None, lon=None, date=None):
    if u is None:
        u = np.zeros(len(z))
    if v is None:
        v = np.zeros(len(z))
    if w is None:
        w = np.zeros(len(z))
    if date is None:
        date = np.array(np.nan)
    if lat is None:
        lat = np.nan
    if lon is None:
        lon = np.nan

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
                    positive='down',
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
        self._horizontal_transform = _roms_horizontal_transform(file_names[0])

    def profile(self, lat, lon, date):
        if isinstance(date, str):
            date = np.datetime64(date)

        file_name, t = self._find_date(date)
        eta, xi = self.horizontal_transform(lat, lon)
        j, i = np.round([eta, xi]).astype(int)

        with _roms_open(file_name) as dset:
            ddset = dset.isel(
                ocean_time=t,
                eta_rho=j,
                xi_rho=i,
                eta_u=j,
                xi_u=range(i-1, i+1),
                eta_v=range(j-1, j+1),
                xi_v=i,
            )

            u = np.flip(ddset.u.mean(dim='xi_u').values)
            v = np.flip(ddset.v.mean(dim='eta_v').values)
            w = np.zeros_like(u)
            salt = np.flip(ddset.salt.values)
            temp = np.flip(ddset.temp.values)
            z = -np.flip(_roms_get_zrho(ddset))

        dens = _seawater_density(temp, salt, z)

        return _create_profile(z, dens, u, v, w, lat, lon, date)

    def horizontal_transform(self, lat, lon):
        return self._horizontal_transform(lat, lon)

    def _find_date(self, date):
        row_idx = np.searchsorted(self.filetab.time_val.values, date)
        file_name = self.filetab.iloc[row_idx].file_name
        time_idx = self.filetab.iloc[row_idx].time_index
        return file_name, time_idx


@contextlib.contextmanager
def _roms_open(file_name, *args, **kwargs):
    try:
        with xr.open_dataset(file_name, *args, **kwargs) as dset:
            yield dset
    except ValueError:
        yield file_name


def _seawater_density(temp, salt, depth):
    """Simplified equation of state, from the NEMO documentation"""

    rho_0 = 1026         # seawater reference density
    a0 = 1.6550e-1       # linear thermal expansion coeff.
    b0 = 7.6554e-1       # linear haline expansion coeff.
    lambda1 = 5.9520e-2  # cabbeling coeff. in T^2
    lambda2 = 5.4914e-4  # cabbeling coeff. in S^2
    nu = 2.4341e-3       # cabbeling coeff. in T S
    mu1 = 1.4970e-4      # thermobaric coeff. in T
    mu2 = 1.1090e-5      # thermobaric coeff. in S

    t_a = temp - 10
    s_a = salt - 35
    z = depth

    da_times_rho0 = (
            -a0 * (1 + 0.5 * lambda1 * t_a + mu1 * z) * t_a
            + b0 * (1 - 0.5 * lambda2 * s_a - mu2 * z) * s_a
            - nu * t_a * s_a
    )

    return rho_0 + da_times_rho0


def _roms_get_zrho(dset):
    h = dset.h.values.ravel()
    hc = dset.hc.values
    c = dset.Cs_r.values
    vt = dset.Vtransform.values

    s = -1.0 + (0.5 + np.arange(len(c))) / len(c)
    out_shape = (len(c), ) + dset.h.shape

    if vt == 1:  # Default transform by Song and Haidvogel
        A = hc * (s - c)[:, None]
        B = np.outer(c, h)
        z = A + B

    elif vt == 2:  # New transform by Shchepetkin
        n = hc * s[:, None] + np.outer(c, h)
        d = 1.0 + hc / h
        z = n / d

    else:
        raise ValueError("Unknown Vtransform")

    return xr.DataArray(
        data=z.reshape(out_shape),
        dims=dset.Cs_r.dims + dset.h.dims,
        name="z_rho",
    )


def sdepth(H, Hc, C, stagger="rho", Vtransform=1):
    H = np.asarray(H)
    Hshape = H.shape  # Save the shape of H
    H = H.ravel()  # and make H 1D for easy shape maniplation
    C = np.asarray(C)
    N = len(C)
    outshape = (N,) + Hshape  # Shape of output
    if stagger == "rho":
        S = -1.0 + (0.5 + np.arange(N)) / N  # Unstretched coordinates
    elif stagger == "w":
        S = np.linspace(-1.0, 0.0, N)
    else:
        raise ValueError("stagger must be 'rho' or 'w'")

    if Vtransform == 1:  # Default transform by Song and Haidvogel
        A = Hc * (S - C)[:, None]
        B = np.outer(C, H)
        return (A + B).reshape(outshape)

    if Vtransform == 2:  # New transform by Shchepetkin
        N = Hc * S[:, None] + np.outer(C, H)
        D = 1.0 + Hc / H
        return (N / D).reshape(outshape)

    # else:
    raise ValueError("Unknown Vtransform")


def _roms_horizontal_transform(file_name):
    from .utils import bilin_inv

    with _roms_open(file_name) as dset:
        lat_rho = dset.variables['lat_rho'].values
        lon_rho = dset.variables['lon_rho'].values

    previous_latlon = (None, None)
    previous_yx = [None, None]

    def latlon_to_yx(lat, lon):
        if previous_latlon != (lat, lon):
            y, x = bilin_inv(lat, lon, lat_rho, lon_rho)
            previous_yx[0] = y
            previous_yx[1] = x
        return previous_yx

    return latlon_to_yx


def _filetab_create(file_names):
    if isinstance(file_names, str):
        file_names = [file_names]

    import pandas as pd

    file_index_series = []
    file_name_series = []
    time_index_series = []
    time_val_series = []
    for i, fname in enumerate(file_names):
        with _roms_open(fname) as dset:
            t = dset.variables['ocean_time'].values.astype('datetime64[us]').tolist()
        time_val_series += t
        time_index_series += list(range(len(t)))
        file_name_series += [fname] * len(t)
        file_index_series += [i] * len(t)

    # file_name_series may be an object list instead of string list
    # We convert it to numpy array in a safe manner
    file_name_series_numpy = np.array([None] * len(file_name_series), dtype=object)
    for i, fname in enumerate(file_name_series):
        file_name_series_numpy[i] = fname

    df = pd.DataFrame(
        dict(
            file_name=file_name_series_numpy,
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
