from effluent import roms
import pytest
import numpy as np
from pathlib import Path


FORCING_1 = str(Path(__file__).parent.joinpath('forcing_1.nc'))
FORCING_glob = str(Path(__file__).parent.joinpath('forcing_?.nc'))


class Test_open_dataset:
    def test_fails_if_no_files(self):
        with pytest.raises(ValueError):
            with roms.open_dataset(file='this_is_no_file'):
                pass

    def test_correct_dims_when_multi_file(self):
        with roms.open_dataset(file=FORCING_glob) as dset:
            assert dset.dims['ocean_time'] == 8
            assert dset['zeta'].dims == ('ocean_time', 'eta_rho', 'xi_rho')
            assert dset['h'].dims == ('eta_rho', 'xi_rho')
            assert dset['hc'].dims == ()

    def test_can_add_depth_values(self):
        with roms.open_dataset(file=FORCING_1, z_rho=True) as dset:
            assert dset['z_rho_star'].dims == ('s_rho', 'eta_rho', 'xi_rho')
            assert dset['z_rho'].dims == ('ocean_time', 's_rho', 'eta_rho', 'xi_rho')

            depths = dset['z_rho_star'].values
            assert np.all(depths < 0)
            assert np.all(depths > -dset['h'].values)

    def test_can_add_dens_values(self):
        with roms.open_dataset(file=FORCING_1, dens=True) as dset:
            assert dset['dens'].dims == ('ocean_time', 's_rho', 'eta_rho', 'xi_rho')

            dens = dset['dens'].values
            assert np.nanmin(dens) > 900
            assert np.nanmax(dens) < 1100
