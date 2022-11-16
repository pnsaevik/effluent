import xarray as xr
import netCDF4 as nc
import numpy as np
from effluent import model


class Test_write_xr_to_nc:
    def test_writes_data_var_values(self):
        nc_dset = nc.Dataset(filename='test_writes_data_vars', mode='w', diskless=True)
        xr_dset = xr.Dataset(
            data_vars=dict(
                a=xr.Variable(('y', 'x'), np.arange(12).reshape((4, 3)))
            ),
        )
        model.write_xr_to_nc(xr_dset, nc_dset)
        assert nc_dset.variables['a'][:].tolist() == xr_dset['a'].values.tolist()
