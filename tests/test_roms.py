# Extra top-level import to address bug in netCDF4 library
# noinspection PyUnresolvedReferences
import netCDF4

from effluent import roms
import numpy as np
from pathlib import Path
import xarray as xr
import pytest


FIXTURES_DIR = Path(__file__).parent.joinpath('fixtures')
FORCING_1 = str(FIXTURES_DIR / 'forcing_1.nc')
FORCING_glob = str(FIXTURES_DIR / 'forcing_?.nc')


@pytest.fixture(scope='module')
def dset1():
    with xr.open_dataset(FORCING_1) as dset:
        yield dset


class Test_select_xy:
    def test_returns_single_point(self, dset1):
        ds = roms.select_xy(dset1, x=2, y=3)
        assert ds.temp.dims == ('ocean_time', 's_rho')

    def test_returns_correct_values(self, dset1):
        ds = roms.select_xy(dset1, x=2, y=3)
        assert 10.608 < float(ds.temp[0, 0]) < 10.610
        assert 0.018 < float(ds.u[0, 0]) < 0.020
        assert -0.044 < float(ds.v[0, 0]) < -0.042


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


class Test_load_location:
    def test_correct_profile_data_sorted_by_increasing_depth(self):
        dset = roms.load_location(FORCING_glob, lat=59.03, lon=5.68, az=0)
        varnames = ['time', 'depth', 'u', 'v', 'salt', 'temp', 'dens']
        try:
            df = dset[varnames].to_dataframe()

        finally:
            dset.close()

        df = df.reset_index()
        df = df[varnames]
        txt = df.to_csv(index=False, float_format="%.02f", lineterminator='\n')

        fname = FIXTURES_DIR / 'profile.csv'
        with open(fname, encoding='utf-8') as fp:
            expected = fp.read()

        assert txt == expected
