"""
The module contains functions for working with ROMS datasets
"""

from __future__ import annotations

import numpy as np
import glob
import xarray as xr
import effluent.eos
import effluent.numerics
import logging


logger = logging.getLogger(__name__)


def open_location(file, lat, lon, az) -> xr.Dataset:
    """
    Open ROMS dataset at specific location

    The output coordinates are 'time' and 'depth'. Fields are interpolated to the
    desired position, and the depth of each vertical level is computed. Current
    velocities are rotated according to the specified azimuthal orientation.

    :param file: Name of ROMS file(s), or wildcard pattern
    :param lat: Latitude of location
    :param lon: Longitude of location
    :param az: Azimuthal orientation of u velocity (0 is north, 90 is east)
    :return: An xarray.Dataset object
    """

    # Find all files
    if isinstance(file, str):
        fnames = sorted(glob.glob(file))
    else:
        fnames = file

    if len(fnames) == 0:
        raise ValueError(f'No files found: "{fnames}"')

    # Compute depth info
    logger.info(f'Compute depths from {fnames[0]}')
    with xr.open_dataset(fnames[0]) as dset:
        ddset = dset[['Vtransform', 'hc', 'h', 'Cs_r', 'lat_rho', 'lon_rho', 'u', 'v']].isel(ocean_time=0)
        dset_point = interpolate_latlon(ddset, lat, lon)
        zrho_star = compute_zrho_star(dset_point)

    # Extract profile info for each dataset
    profile_dsets = []
    for fname in fnames:
        logger.info(f'Open file {fname}')
        with xr.open_dataset(fname) as dset:
            logger.info(f'Horizontal interpolation')
            dset = dset[['u', 'v', 'temp', 'salt', 'angle']]
            dset = interpolate_xy(dset, dset_point.xi_rho.values, dset_point.eta_rho.values)

            logger.info(f'Rotate velocity vectors, compute density')
            u = compute_azimuthal_vel(dset, az * (np.pi / 180))
            v = compute_azimuthal_vel(dset, (az + 90) * (np.pi / 180))
            dset = dset.assign(u=u, v=v)

            dset = dset.assign(z_rho_star=zrho_star)
            dset = dset.assign(dens=compute_dens(dset))
            dset = dset.rename(z_rho_star='depth', ocean_time='time')
            dset = dset.assign_coords(depth=-dset['depth'])
            dset = dset.swap_dims({'s_rho': 'depth'})

            profile_dsets.append(dset)

            logger.debug(f'Close file {fname}')

    # Concatenate datasets
    logger.info('Concatenate datasets')
    dset_combined = xr.concat(
        objs=profile_dsets,
        dim='time',
        data_vars='minimal',
        coords='minimal',
        compat='override',
        combine_attrs='override',
    )

    return dset_combined


def compute_zrho(dset: xr.Dataset) -> xr.Dataset:
    """
    Compute z_rho variable from a ROMS dataset

    The z_rho variable is negative depth, with tidal variation

    :param dset: ROMS dataset
    :return: The z_rho variable
    """
    vtrans = dset['Vtransform']
    z_rho_star = dset['z_rho_star']

    if vtrans == 1:
        z_rho = z_rho_star + dset.zeta * (1 + z_rho_star / dset.h)
    elif vtrans == 2:
        z_rho = dset.zeta + z_rho_star * (dset.zeta / dset.h + 1)
    else:
        raise ValueError(f'Unknown Vtransform: {vtrans}')

    dims = [d for d in ['ocean_time', 's_rho', 'eta_rho', 'xi_rho'] if d in z_rho.dims]
    z_rho = z_rho.transpose(*dims)
    z_rho.name = 'z_rho'
    return z_rho


def compute_zrho_star(dset: xr.Dataset) -> xr.DataArray:
    """
    Compute z_rho_star variable from a ROMS dataset

    The z_rho_star variable is negative depth, without tidal variation

    :param dset: ROMS dataset
    :return: The z_rho_star variable
    """
    vtrans = dset['Vtransform']

    if vtrans == 1:
        z_rho_star = dset.hc * (dset.s_rho - dset.Cs_r) + dset.Cs_r * dset.h
    elif vtrans == 2:
        z_rho_0 = (dset.hc * dset.s_rho + dset.Cs_r * dset.h) / (dset.hc + dset.h)
        z_rho_star = z_rho_0 * dset.h
    else:
        raise ValueError(f'Unknown Vtransform: {vtrans}')

    dims = [d for d in ['s_rho', 'eta_rho', 'xi_rho'] if d in z_rho_star.dims]
    z_rho_star = z_rho_star.transpose(*dims)
    z_rho_star.name = 'z_rho_star'
    return z_rho_star


def compute_dens(dset: xr.Dataset) -> xr.DataArray:
    """
    Compute variable ``dens`` from a ROMS dataset

    :param dset: ROMS dataset
    :return: New dataset with ``dens`` added
    """
    dens = effluent.eos.roms_rho(dset.temp, dset.salt, dset.z_rho_star)
    dens.name = 'dens'
    return dens


def interpolate_latlon(dset: xr.Dataset, lat, lon) -> xr.Dataset:
    """
    Interpolate fields in ROMS dataset

    The function uses bilinear interpolation for regular field variables, and
    unidirectional interpolation (which preserves divergence) for the ``u`` and ``v``
    variables.

    :param dset: ROMS dataset
    :param lat: The latitude
    :param lon: The longitude
    :return: New dataset with all variables interpolated to the specified location
    """
    lat_rho = dset.lat_rho.values
    lon_rho = dset.lon_rho.values

    y, x = effluent.numerics.bilin_inv(lat, lon, lat_rho, lon_rho)

    return interpolate_xy(dset, x, y)


def interpolate_xy(dset: xr.Dataset, x, y) -> xr.Dataset:
    """
    Interpolate fields in ROMS dataset

    The function uses bilinear interpolation for regular field variables, and
    unidirectional interpolation (which preserves divergence) for the ``u`` and ``v``
    variables.

    :param dset: ROMS dataset
    :param x: The dataset x coordinate
    :param y: The dataset y coordinate
    :return: New dataset with all variables interpolated to the specified location
    """
    x_min = 0.5
    y_min = 0.5
    x_max = dset.dims['xi_rho'] - 1.5
    y_max = dset.dims['eta_rho'] - 1.5
    x = np.clip(x, x_min, x_max)
    y = np.clip(y, y_min, y_max)

    cvars = {'xi_rho', 'xi_u', 'xi_v', 'eta_rho', 'eta_u', 'eta_v'}
    dset = dset.drop_vars(cvars.intersection(dset.variables))
    dset = dset.interp(
        xi_rho=x,
        eta_rho=y,
        xi_u=x - 0.5,
        eta_u=int(y + 0.5),
        xi_v=int(x + 0.5),
        eta_v=y - 0.5,
    )

    return dset


def compute_azimuthal_vel(dset: xr.Dataset, az) -> xr.DataArray:
    """
    Compute directional current velocity

    :param dset: ROMS dataset
    :param az: The direction in which to measure the current
    :return: An xarray.DataArray representing the current velocity
    """

    assert dset.angle.units == "radians"

    u = dset.u
    v = dset.v
    theta = az + np.pi / 2 - dset.angle
    return u * np.cos(theta) + v * np.sin(theta)
