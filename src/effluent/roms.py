import xarray as xr
import numpy as np


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
