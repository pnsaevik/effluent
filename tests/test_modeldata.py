import pandas as pd
import pytest

from effluent import modeldata
import numpy as np
from pathlib import Path
import xarray as xr


root = Path(__file__).parent
file_name_1 = str(root / 'forcing_1.nc')
file_name_2 = str(root / 'forcing_2.nc')


@pytest.fixture(scope='module')
def dataset_1():
    return xr.load_dataset(file_name_1)


@pytest.fixture(scope='module')
def dataset_2():
    return xr.load_dataset(file_name_2)


class Test_from_profile:
    def test_returns_dataset_with_standardized_names(self):
        m = modeldata.Modeldata.from_profile(
            z=[0, 10, 50],
            u=[.3, .2, .1],
            v=[.03, .02, .01],
            dens=[1023, 1024, 1025],
        )

        p = m.profile()

        assert set(p.variables) == {'z', 'dens', 'u', 'v', 'w', 'lat', 'lon', 'date'}


class Test_from_roms:
    def test_returns_dataset_with_standardized_names(self, dataset_1, dataset_2):
        m = modeldata.Modeldata.from_roms([dataset_1, dataset_2])
        p = m.profile(59.02619, 5.67831, '2015-09-07T03')

        assert set(p.variables) == {'z', 'dens', 'u', 'v', 'w', 'lat', 'lon', 'date'}

    def test_returned_profile_has_positive_and_increasing_depth(self, dataset_1, dataset_2):
        m = modeldata.Modeldata.from_roms([dataset_1, dataset_2])
        p = m.profile(59.02619, 5.67831, '2015-09-07T03')

        assert all(d > 0 for d in p.z.values)
        assert all(d > 0 for d in np.diff(p.z.values))


class Test_filetab_create:
    def test_returns_merged_table_of_times(self):
        df = modeldata._filetab_create([file_name_1, file_name_2])
        assert df.file_index.tolist() == [0, 0, 0, 0, 1, 1, 1, 1]
        assert df.file_name.tolist() == [file_name_1] * 4 + [file_name_2] * 4
        assert df.time_index.tolist() == [0, 1, 2, 3, 0, 1, 2, 3]
        assert isinstance(df.time_val[0], pd.Timestamp)
        assert df.time_val.values.astype('datetime64[h]').astype(str).tolist() == [
            '2015-09-07T01', '2015-09-07T02', '2015-09-07T03', '2015-09-07T04',
            '2015-09-07T05', '2015-09-07T06', '2015-09-07T07', '2015-09-07T08',
        ]

    def test_returns_merged_table_of_times_when_input_is_object(self, dataset_1, dataset_2):
        df = modeldata._filetab_create([dataset_1, dataset_2])
        assert df.file_index.tolist() == [0, 0, 0, 0, 1, 1, 1, 1]
        assert df.file_name.tolist() == [dataset_1] * 4 + [dataset_2] * 4
        assert df.time_index.tolist() == [0, 1, 2, 3, 0, 1, 2, 3]
        assert df.time_val.values.astype('datetime64[h]').astype(str).tolist() == [
            '2015-09-07T01', '2015-09-07T02', '2015-09-07T03', '2015-09-07T04',
            '2015-09-07T05', '2015-09-07T06', '2015-09-07T07', '2015-09-07T08',
        ]


class Test_horizontal_transform:
    def test_returns_transform_from_latlon_to_ji(self, dataset_1):
        horz_trans = modeldata._roms_horizontal_transform(dataset_1)
        j, i = horz_trans(59.02619, 5.67831)
        assert i.round(2) == 3.00
        assert j.round(2) == 2.00


class Test_from_date_to_index:
    def test_returns_file_name_and_time_index(self, dataset_1, dataset_2):
        m = modeldata.ModeldataRoms([dataset_1, dataset_2])
        fname, tidx = m._find_date(np.datetime64('2015-09-07T07'))

        assert fname is dataset_2
        assert tidx == 2


class Test_roms_get_zrho:
    def test_returns_z_array(self, dataset_1):
        zrho = modeldata._roms_get_zrho(dataset_1)
        assert isinstance(zrho, xr.DataArray)


