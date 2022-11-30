from effluent import roms
import xarray as xr
import pytest


class Test_open_dataset:
    def test_fails_if_no_files(self):
        with pytest.raises(ValueError):
            with roms.open_dataset(file='this_is_no_file'):
                pass

    def test_returns_dataset_if_single_file(self):
        with roms.open_dataset(file='forcing_1.nc') as dset:
            assert isinstance(dset, xr.Dataset)

    def test_correct_dims_when_multi_file(self):
        with roms.open_dataset(file='forcing_?.nc') as dset:
            assert dset.dims['ocean_time'] == 8
            assert dset['zeta'].dims == ('ocean_time', 'eta_rho', 'xi_rho')
            assert dset['h'].dims == ('eta_rho', 'xi_rho')
            assert dset['hc'].dims == ()
