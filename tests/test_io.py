import io
import uuid

import netCDF4 as nc
import numpy as np
import pytest
import xarray as xr

import effluent.io


class Test_write_xr_to_nc:
    @pytest.fixture()
    def nc_dset(self):
        fname = uuid.uuid4()
        with nc.Dataset(filename=fname, mode='w', diskless=True) as dset:
            yield dset

    def test_writes_data_var_values(self, nc_dset):
        xr_dset = xr.Dataset(
            data_vars=dict(
                a=xr.Variable(('y', 'x'), np.arange(12).reshape((4, 3)))
            ),
        )
        effluent.io.write_xr_to_nc(xr_dset, nc_dset)
        assert nc_dset.variables['a'][:].tolist() == xr_dset['a'].values.tolist()

    def test_writes_data_coord_values(self, nc_dset):
        xr_dset = xr.Dataset(coords=dict(a=xr.Variable('a', np.arange(12))))
        effluent.io.write_xr_to_nc(xr_dset, nc_dset)
        assert nc_dset.variables['a'][:].tolist() == xr_dset['a'].values.tolist()

    def test_writes_variable_attrs(self, nc_dset):
        xr_dset = xr.Dataset(
            data_vars=dict(
                a=xr.Variable('x', np.arange(5), attrs=dict(units='m'))
            )
        )
        effluent.io.write_xr_to_nc(xr_dset, nc_dset)
        assert nc_dset.variables['a'].units == xr_dset['a'].attrs['units']

    def test_writes_dataset_attrs(self, nc_dset):
        xr_dset = xr.Dataset(
            data_vars=dict(a=xr.Variable('x', np.arange(5))),
            attrs=dict(date='2000-01-01'),
        )
        effluent.io.write_xr_to_nc(xr_dset, nc_dset)
        assert nc_dset.date == xr_dset.attrs['date']

    def test_writes_unlimited_dims(self, nc_dset):
        xr_dset = xr.Dataset(
            data_vars=dict(a=xr.Variable('x', np.arange(5))),
        )
        xr_dset.encoding['unlimited_dims'] = ['x']
        effluent.io.write_xr_to_nc(xr_dset, nc_dset)
        assert nc_dset.dimensions['x'].isunlimited()

    def test_writes_datetimes_as_seconds_since_epoch(self, nc_dset):
        datetimes = np.array(['1970-01-01', '1970-01-01T01']).astype('datetime64')
        xr_dset = xr.Dataset(data_vars=dict(t=xr.Variable('t', datetimes)))

        effluent.io.write_xr_to_nc(xr_dset, nc_dset)
        assert nc_dset.variables['t'][:].tolist() == [0, 3600]
        assert nc_dset.variables['t'].units == 'seconds since 1970-01-01'
        assert nc_dset.variables['t'].dtype == np.dtype('i8')


class Test_append_xr_to_nc:
    @pytest.fixture()
    def nc_dset(self):
        fname = uuid.uuid4()
        with nc.Dataset(filename=fname, mode='w', diskless=True) as dset:
            dset.createDimension('x', 4)
            dset.createDimension('y', None)
            dset.createVariable('a', 'i4', ('y', 'x'))
            dset.createVariable('b', 'i4', 'x')
            dset.createVariable('t', 'i8', 'y')
            dset.variables['t'].units = 'hours since 1970-01-01'
            dset.variables['t'].calendar = 'proleptic_gregorian'
            dset.variables['a'][:2, :] = 0
            dset.variables['b'][:] = 1
            dset.variables['t'][:2] = -1
            yield dset

    def test_appends_variable_data_to_dataset_if_unlim_dims(self, nc_dset):
        xr_dset = xr.Dataset(
            data_vars=dict(
                a=xr.Variable(('y', 'x'), np.arange(12).reshape((3, 4)))
            ),
        )
        effluent.io.append_xr_to_nc(xr_dset, nc_dset)
        assert nc_dset.variables['a'][:].tolist() == [
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 1, 2, 3],
            [4, 5, 6, 7],
            [8, 9, 10, 11],
        ]

    def test_ignores_variables_without_unlim_dims(self, nc_dset):
        xr_dset = xr.Dataset(data_vars=dict(b=xr.Variable('x', np.arange(3))))
        assert nc_dset['b'][:].tolist() == [1, 1, 1, 1]
        effluent.io.append_xr_to_nc(xr_dset, nc_dset)
        assert nc_dset['b'][:].tolist() == [1, 1, 1, 1]

    def test_writes_datetimes_using_preexisting_units(self, nc_dset):
        datetimes = np.array(['1970-01-01', '1970-01-01T05']).astype('datetime64')
        xr_dset = xr.Dataset(data_vars=dict(t=xr.Variable('y', datetimes)))
        effluent.io.append_xr_to_nc(xr_dset, nc_dset)
        assert nc_dset.variables['t'][:].tolist() == [-1, -1, 0, 5]


class Test_Pipe_from_config:
    def test_explicit_mapping_with_arrays(self):
        posix = np.array([0, 1200])
        time = np.datetime64('1970-01-01') + posix.astype('timedelta64[s]')
        conf = dict(
            time=time.astype(object), flow=[2, 3], dens=[4, 5], diam=[6, 7],
            depth=[8, 9], decline=[10, 11],
        )
        p = effluent.io.Pipe.from_config(conf)
        dset = p.select(time=np.datetime64('1970-01-01 00:10'))
        assert dset.dens.values == 4.5

    def test_explicit_mapping_with_singlenums(self):
        posix = np.array([0, 1200])
        time = np.datetime64('1970-01-01') + posix.astype('timedelta64[s]')
        conf = dict(
            time=time.astype(object), flow=[2, 3], dens=4, diam=6, depth=[8, 9],
            decline=10,
        )
        p = effluent.io.Pipe.from_config(conf)
        dset = p.select(time=np.datetime64('1970-01-01 00:10'))
        assert dset.dens.values == 4

    def test_can_convert_salt_and_temp_to_dens(self):
        posix = np.array([0, 1200])
        time = np.datetime64('1970-01-01') + posix.astype('timedelta64[s]')
        conf = dict(
            time=time.astype(object), flow=[2, 3], diam=[6, 7],
            depth=[8, 9], decline=[10, 11], salt=[30, 32], temp=[4, 5],
        )
        p = effluent.io.Pipe.from_config(conf)
        dset = p.select(time=np.datetime64('1970-01-01 00:10'))
        assert 1024 < dset.dens.values.item() < 1025

    def test_can_have_single_time_entry(self):
        posix = np.array([0])
        time = np.datetime64('1970-01-01') + posix.astype('timedelta64[s]')
        conf = dict(
            time=time.astype(object), flow=[2], dens=[4], diam=[6], depth=[8],
            decline=[10],
        )
        p = effluent.io.Pipe.from_config(conf)
        dset = p.select(time=np.datetime64('1970-01-01 00:10'))
        assert dset.dens.values.item() == 4

    def test_csv_file(self):
        buf = io.StringIO("""
                        time, flow, dens, diam, depth, decline
            1970-01-01 00:00,    1,    3,    5,     7,       9
            1970-01-01 00:20,    2,    4,    6,     8,      10
        """)
        conf = dict(csv=dict(file=buf))
        p = effluent.io.Pipe.from_config(conf)
        dset = p.select(time=np.datetime64('1970-01-01 00:10'))
        assert dset.dens.values == 3.5

    def test_nc_file(self):
        posix = np.array([0, 1200])
        time = np.datetime64('1970-01-01') + posix.astype('timedelta64[s]')
        buf = xr.Dataset(
            coords=dict(time=time),
            data_vars=dict(
                flow=xr.Variable('time', [1, 2]),
                dens=xr.Variable('time', [3, 4]),
                diam=xr.Variable('time', [5, 6]),
                depth=xr.Variable('time', [7, 8]),
                decline=xr.Variable('time', [9, 10]),
            ),
        ).to_netcdf()
        conf = dict(nc=dict(file=buf))
        p = effluent.io.Pipe.from_config(conf)
        dset = p.select(time=np.datetime64('1970-01-01 00:10'))
        assert dset.dens.values == 3.5


class Test_Ambient_from_config:
    def test_explicit_mapping_with_arrays(self):
        posix = np.array([0, 1200])
        time = np.datetime64('1970-01-01') + posix.astype('timedelta64[s]')
        conf = dict(
            time=time.astype(object),
            depth=[0, 10, 20],
            coflow=[[0, 1, 2], [3, 4, 5]],
            crossflow=[[6, 7, 8], [9, 0, 1]],
            dens=[[2, 3, 4], [5, 6, 7]],
        )
        p = effluent.io.Ambient.from_config(conf)
        dset = p.select(time=np.datetime64('1970-01-01 00:10'))
        assert dset.dens.values.tolist() == [3.5, 4.5, 5.5]

    def test_can_have_single_time_entry(self):
        posix = np.array([0])
        time = np.datetime64('1970-01-01') + posix.astype('timedelta64[s]')
        conf = dict(
            time=time.astype(object),
            depth=[0, 10, 20],
            coflow=[[0, 1, 2]],
            crossflow=[[6, 7, 8]],
            dens=[[2, 3, 4]],
        )
        p = effluent.io.Ambient.from_config(conf)
        dset = p.select(time=np.datetime64('1970-01-01 00:10'))
        assert dset.dens.values.tolist() == [2, 3, 4]

    def test_explicit_mapping_with_singlenums(self):
        posix = np.array([0, 1200])
        time = np.datetime64('1970-01-01') + posix.astype('timedelta64[s]')
        conf = dict(
            time=time.astype(object),
            depth=[0, 10, 20],
            coflow=[[0, 1, 2], [3, 4, 5]],
            crossflow=[[6, 7, 8], [9, 0, 1]],
            dens=3,
        )
        p = effluent.io.Ambient.from_config(conf)
        dset = p.select(time=np.datetime64('1970-01-01 00:10'))
        assert dset.dens.values.tolist() == [3, 3, 3]

    def test_can_convert_salt_and_temp_to_dens(self):
        posix = np.array([0, 1200])
        time = np.datetime64('1970-01-01') + posix.astype('timedelta64[s]')
        conf = dict(
            time=time.astype(object),
            depth=[0, 10, 20],
            coflow=[[0, 1, 2], [3, 4, 5]],
            crossflow=[[6, 7, 8], [9, 0, 1]],
            temp=[[2, 3, 4], [5, 6, 7]],
            salt=[[22, 23, 24], [25, 26, 27]],
        )
        p = effluent.io.Ambient.from_config(conf)
        dset = p.select(time=np.datetime64('1970-01-01 00:10'))
        assert 1019 < dset.dens.values.mean() < 1020

    def test_csv_file(self):
        buf = io.StringIO("""
                        time, depth, coflow, crossflow, dens
            1970-01-01 00:00,     0,      0,         6,    2
            1970-01-01 00:00,    10,      1,         7,    3
            1970-01-01 00:00,    20,      2,         8,    4
            1970-01-01 00:20,     0,      3,         9,    5
            1970-01-01 00:20,    10,      4,         0,    6
            1970-01-01 00:20,    20,      5,         1,    7
        """)
        conf = dict(csv=dict(file=buf))
        p = effluent.io.Ambient.from_config(conf)
        dset = p.select(time=np.datetime64('1970-01-01 00:10'))
        assert dset.dens.values.tolist() == [3.5, 4.5, 5.5]

    def test_nc_file(self):
        posix = np.array([0, 1200])
        time = np.datetime64('1970-01-01') + posix.astype('timedelta64[s]')
        buf = xr.Dataset(
            coords=dict(
                time=time,
                depth=[0, 10, 20],
            ),
            data_vars=dict(
                coflow=xr.Variable(('time', 'depth'), [[0, 1, 2], [3, 4, 5]]),
                crossflow=xr.Variable(('time', 'depth'), [[6, 7, 8], [9, 0, 1]]),
                dens=xr.Variable(('time', 'depth'), [[2, 3, 4], [5, 6, 7]]),
            )
        ).to_netcdf()
        conf = dict(nc=dict(file=buf))
        p = effluent.io.Ambient.from_config(conf)
        dset = p.select(time=np.datetime64('1970-01-01 00:10'))
        assert dset.dens.values.tolist() == [3.5, 4.5, 5.5]

    def test_roms_file(self):
        from pathlib import Path
        FORCING_glob = str(Path(__file__).parent.joinpath('fixtures/forcing_?.nc'))

        conf = dict(roms=dict(file=FORCING_glob, latitude=59.03, longitude=5.68, azimuth=0))
        with effluent.io.Ambient.from_config(conf) as p:
            dset = p.select(time=np.datetime64('2015-09-07T07'))
            assert len(dset.depth) > 0
            assert dset.dens.dims == ('depth', )


class Test_Output_from_config:
    @pytest.fixture()
    def result(self):
        return xr.Dataset(
            data_vars=dict(
                x=xr.Variable('t', [1, 2]),
                y=xr.Variable('t', [3, 4]),
                z=xr.Variable('t', [5, 6]),
                u=xr.Variable('t', [7, 8]),
                v=xr.Variable('t', [9, 1]),
                w=xr.Variable('t', [2, 3]),
                density=xr.Variable('t', [4, 5]),
                radius=xr.Variable('t', [6, 7]),
            ),
            coords=dict(
                t=xr.Variable('t', [1000, 2000]),
            ),
        )

    def test_option_csv_file(self, result):
        buf = io.StringIO()
        conf = dict(
            csv=dict(file=buf),
        )

        with effluent.io.Output.from_config(conf) as out:
            out.write(time=0, result=result)
            out.write(time=1, result=result)
            txt = buf.getvalue()

        assert txt.replace('\r', '') == (
            "release_time,t,x,y,z,u,v,w,density,radius\n"
            "0,1000,1,3,5,7,9,2,4,6\n"
            "0,2000,2,4,6,8,1,3,5,7\n"
            "1,1000,1,3,5,7,9,2,4,6\n"
            "1,2000,2,4,6,8,1,3,5,7\n"
        )

    def test_option_nc_file(self, result):
        buf = xr.Dataset()
        conf = dict(
            nc=dict(file=buf),
        )

        with effluent.io.Output.from_config(conf) as out:
            out.write(time=0, result=result)
            out.write(time=1, result=result)

        assert buf.z.positive == 'down'
        assert buf.z.dims == ('release_time', 't')
        assert buf.z.shape == (2, 2)
        assert buf.z.values.tolist() == [[5, 6], [5, 6]]
        assert set(buf.data_vars) == {'x', 'y', 'z', 'u', 'v', 'w', 'density', 'radius'}
        assert set(buf.coords) == {'release_time', 't'}

    def test_option_variables_when_nc_file(self, result):
        buf = xr.Dataset()
        conf = dict(
            variables=['release_time', 'x', 'y', 'z'],
            nc=dict(file=buf),
        )

        with effluent.io.Output.from_config(conf) as out:
            out.write(time=0, result=result)

        assert set(buf.data_vars) == {'x', 'y', 'z'}
        assert set(buf.coords) == {'release_time'}
        assert set(buf.dims) == {'release_time', 't'}

    def test_option_variables_when_csv_file(self, result):
        buf = io.StringIO()
        conf = dict(
            variables=['release_time', 'x', 'z', 'y'],
            csv=dict(file=buf),
        )

        with effluent.io.Output.from_config(conf) as out:
            out.write(time=0, result=result)

            buf.seek(0)
            first_line = buf.readline()

        assert first_line.replace('\r', '') == "release_time,x,z,y\n"
