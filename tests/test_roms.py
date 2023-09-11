from effluent import roms
import pytest
import numpy as np
from pathlib import Path
import xarray as xr


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


class Test_select_latlon:
    def test_returns_single_point(self):
        with roms.open_dataset(file=FORCING_1, z_rho=True) as dset:
            dset2 = roms.interpolate_latlon(dset, lat=59.03, lon=5.68)
            assert dset2.temp.dims == ('ocean_time', 's_rho')


class Test_compute_azimuthal_vel:
    def test_returns_correct_values(self):
        angle_xi = np.pi/3
        angle_vel = np.array(
            [0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi/6, 5*np.pi/6, 8*np.pi/6, 11*np.pi/6]
        )

        dset = xr.Dataset(
            data_vars=dict(
                u=xr.Variable('t', [30] * 8),
                v=xr.Variable('t', [40] * 8),
                angle=xr.Variable('t', [angle_xi] * 8, attrs={
                    'units': 'radians', 'long_name': 'angle between XI axis and EAST'
                }),
            )
        )
        vel = roms.compute_azimuthal_vel(dset, az=angle_vel)
        assert vel.values.round().astype('i4').tolist() == [46, 20, -46, -20, 40, -30, -40, 30]


class Test_open_location:
    def test_correct_vars_and_dims(self):
        dset = xr.Dataset()
        try:
            dset = roms.open_location(FORCING_glob, lat=59.03, lon=5.68, az=0)
            assert dset.dens.dims == ('time', 'depth')
            assert dset.u.dims == ('time', 'depth')
            assert dset.v.dims == ('time', 'depth')

        finally:
            dset.close()
