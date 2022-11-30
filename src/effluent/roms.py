import xarray as xr
import numpy as np
import glob


def open_location(**kwargs):
    return xr.Dataset(
        data_vars=dict(
            time=np.array(['1970-01-01', '1970-01-01T00:20']).astype('datetime64[s]'),
            depth=[0, 10, 20],
        ),
        coords=dict(
            u=xr.Variable(('time', 'depth'), [[0, 1, 2], [3, 4, 5]]),
            v=xr.Variable(('time', 'depth'), [[6, 7, 8], [9, 0, 1]]),
            dens=xr.Variable(('time', 'depth'), [[2, 3, 4], [5, 6, 7]]),
        ),
    )


def open_dataset(file, z_rho=False):
    """
    Open ROMS dataset

    Variables are lazily loaded or computed.

    :param file: Name of ROMS file(s), or wildcard pattern
    :param z_rho: True if rho depths should be added (default: False)
    :return: An xarray.Dataset object
    """
    fnames = sorted(glob.glob(file))
    if len(fnames) == 0:
        raise ValueError(f'No files found: "{fnames}"')

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

    if z_rho:
        dset = add_zrho(dset)

    return dset


def add_zrho(dset):
    vtrans = dset.get('Vtransform', None)
    if vtrans in [1, 2]:
        if dset.Vtransform == 1:
            z_rho_star = dset.hc * (dset.s_rho - dset.Cs_r) + dset.Cs_r * dset.h
            z_rho = z_rho_star + dset.zeta * (1 + z_rho_star / dset.h)
        else:
            z_rho_0 = (dset.hc * dset.s_rho + dset.Cs_r * dset.h) / (dset.hc + dset.h)
            z_rho_star = z_rho_0 * dset.h
            z_rho = dset.zeta + z_rho_0 * (dset.zeta + dset.h)

        return dset.assign_coords(
            z_rho=z_rho.transpose('ocean_time', 's_rho', 'eta_rho', 'xi_rho'),
            z_rho_star=z_rho_star.transpose('s_rho', 'eta_rho', 'xi_rho'),
        )

    else:
        return dset


def add_dens(dset):
    pass
