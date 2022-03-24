from effluent import modeldata
import numpy as np


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


class Test_filetab_create:
    def test_returns_merged_table_of_times(self):
        from pathlib import Path
        root = Path(__file__).parent
        file_names = [str(root / f) for f in ['forcing_1.nc', 'forcing_2.nc']]
        df = modeldata._filetab_create(file_names)
        assert df.file_index.tolist() == [0, 0, 0, 0, 1, 1, 1, 1]
        assert df.file_name.tolist() == [file_names[0]] * 4 + [file_names[1]] * 4
        assert df.time_index.tolist() == [0, 1, 2, 3, 0, 1, 2, 3]
        assert df.time_val.values.astype('datetime64[h]').astype(str).tolist() == [
            '2015-09-07T01', '2015-09-07T02', '2015-09-07T03', '2015-09-07T04',
            '2015-09-07T05', '2015-09-07T06', '2015-09-07T07', '2015-09-07T08',
        ]


class Test_horizontal_transform:
    def test_returns_transform_from_latlon_to_ji(self):
        from pathlib import Path
        root = Path(__file__).parent
        file_name = str(root / 'forcing_1.nc')
        horz_trans = modeldata._roms_horizontal_transform(file_name)
        j, i = horz_trans(59.02619, 5.67831)
        assert i.round(2) == 3.00
        assert j.round(2) == 2.00


class Test_from_date_to_index:
    def test_returns_file_and_time_index(self):
        from pathlib import Path
        root = Path(__file__).parent
        file_names = [str(root / f) for f in ['forcing_1.nc', 'forcing_2.nc']]

        m = modeldata.ModeldataRoms(file_names)
        fidx, tidx = m._from_date_to_index(np.datetime64('2015-09-07T07'))

        assert fidx == 1
        assert tidx == 2
