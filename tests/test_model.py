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

    def test_writes_dataset_attrs(self):
        nc_dset = nc.Dataset(filename='writes_dset_attrs', mode='w', diskless=True)
        xr_dset = xr.Dataset(
            data_vars=dict(a=xr.Variable('x', np.arange(5))),
            attrs=dict(date='2000-01-01'),
        )
        model.write_xr_to_nc(xr_dset, nc_dset)
        assert nc_dset.date == xr_dset.attrs['date']

    def test_writes_unlimited_dims(self):
        nc_dset = nc.Dataset(filename='writes_unlim_dims', mode='w', diskless=True)
        xr_dset = xr.Dataset(
            data_vars=dict(a=xr.Variable('x', np.arange(5))),
        )
        xr_dset.encoding['unlimited_dims'] = ['x']
        model.write_xr_to_nc(xr_dset, nc_dset)
        assert nc_dset.dimensions['x'].isunlimited()
