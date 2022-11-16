import xarray as xr
import netCDF4 as nc
import numpy as np
from effluent import model


class Test_write_xr_to_nc:
    def test_writes_data_var_values(self):
        nc_dset = nc.Dataset(filename='writes_data_vars', mode='w', diskless=True)
        xr_dset = xr.Dataset(
            data_vars=dict(
                a=xr.Variable(('y', 'x'), np.arange(12).reshape((4, 3)))
            ),
        )
        model.write_xr_to_nc(xr_dset, nc_dset)
        assert nc_dset.variables['a'][:].tolist() == xr_dset['a'].values.tolist()

    def test_writes_data_coord_values(self):
        nc_dset = nc.Dataset(filename='writes_data_coord', mode='w', diskless=True)
        xr_dset = xr.Dataset(coords=dict(a=xr.Variable('a', np.arange(12))))
        model.write_xr_to_nc(xr_dset, nc_dset)
        assert nc_dset.variables['a'][:].tolist() == xr_dset['a'].values.tolist()

    def test_writes_variable_attrs(self):
        nc_dset = nc.Dataset(filename='writes_var_attrs', mode='w', diskless=True)
        xr_dset = xr.Dataset(
            data_vars=dict(
                a=xr.Variable('x', np.arange(5), attrs=dict(units='m'))
            )
        )
        model.write_xr_to_nc(xr_dset, nc_dset)
        assert nc_dset.variables['a'].units == xr_dset['a'].attrs['units']
