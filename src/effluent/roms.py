import xarray as xr
import numpy as np
import glob


def open_location(file, latitude, longitude, azimuth):
    """
    Open a ROMS dataset and return data for a specific location

    The output coordinates are 'time' and 'depth'

    :param file: Name of ROMS file(s), or wildcard pattern
    :param latitude: Latitude of location
    :param longitude: Longitude of location
    :param azimuth: Azimuthal orientation of u velocity (0 is north, 90 is east)
    :return: An xarray.Dataset object
    """

    # Select position
    dset = open_dataset(file, z_rho=True, dens=True)
    dset = select_latlon(dset, latitude, longitude)

    # Set coordinates
    dset = dset.rename(z_rho_star='depth', ocean_time='time')
    dset = dset.swap_dims({'s_rho': 'depth'})

    # Rotate velocity
    u = compute_azimuthal_velocity(dset, azimuth * (np.pi / 180))
    v = compute_azimuthal_velocity(dset, (azimuth + 90) * (np.pi / 180))
    dset = dset.assign(u=u, v=v)

    return dset


def open_dataset(file, z_rho=False, dens=False):
    """
    Open ROMS dataset

    Variables are lazily loaded or computed.

    :param file: Name of ROMS file(s), or wildcard pattern
    :param z_rho: True if rho depths should be added (default: False)
    :param dens: True if density should be added (implies z_rho, default: False)
    :return: An xarray.Dataset object
    """
    fnames = sorted(glob.glob(file))
    if len(fnames) == 0:
        raise ValueError(f'No files found: "{fnames}"')

    if len(fnames) == 1:
        dset = xr.open_dataset(fnames[0])
    else:
        dset = xr.open_mfdataset(
            paths=fnames,
            chunks={'ocean_time': 1},
            concat_dim='ocean_time',
            compat='override',
            data_vars='minimal',
            coords='minimal',
            combine='nested',
            join='override',
            combine_attrs='override',
        )

    if z_rho or dens:
        dset = add_zrho(dset)

    if dens:
        dset = add_dens(dset)

    return dset


def add_zrho(dset):
    vtrans = dset['Vtransform']

    if vtrans == 1:
        z_rho_star = dset.hc * (dset.s_rho - dset.Cs_r) + dset.Cs_r * dset.h
        z_rho = z_rho_star + dset.zeta * (1 + z_rho_star / dset.h)
    elif vtrans == 2:
        z_rho_0 = (dset.hc * dset.s_rho + dset.Cs_r * dset.h) / (dset.hc + dset.h)
        z_rho_star = z_rho_0 * dset.h
        z_rho = dset.zeta + z_rho_0 * (dset.zeta + dset.h)
    else:
        raise ValueError(f'Unknown Vtransform: {vtrans}')

    return dset.assign_coords(
        z_rho=z_rho.transpose('ocean_time', 's_rho', 'eta_rho', 'xi_rho'),
        z_rho_star=z_rho_star.transpose('s_rho', 'eta_rho', 'xi_rho'),
    )


def add_dens(dset):
    from effluent.eos import roms_rho
    dens = roms_rho(dset.temp, dset.salt, dset.z_rho_star)
    return dset.assign_coords(dens=dens)


def select_latlon(dset, lat, lon):
    from .numerics import bilin_inv

    lat_rho = dset.lat_rho.values
    lon_rho = dset.lon_rho.values

    y, x = bilin_inv(lat, lon, lat_rho, lon_rho)

    x_min = 0.5
    y_min = 0.5
    x_max = dset.dims['xi_rho'] - 1.5
    y_max = dset.dims['eta_rho'] - 1.5
    x = np.clip(x, x_min, x_max)
    y = np.clip(y, y_min, y_max)

    dset = dset.interp(
        xi_rho=x,
        eta_rho=y,
        xi_u=x - 0.5,
        eta_u=int(y + 0.5),
        xi_v=int(x + 0.5),
        eta_v=y - 0.5,
    )

    return dset


def compute_azimuthal_velocity(dset, azimuth):
    assert dset.angle.units == "radians"

    u = dset.u
    v = dset.v
    theta = azimuth + np.pi/2 - dset.angle
    return u * np.cos(theta) + v * np.sin(theta)
