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


def open_dataset(file):
    fnames = sorted(glob.glob(file))
    if len(fnames) == 0:
        raise ValueError(f'No files found: "{fnames}"')

    return xr.open_mfdataset(
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
