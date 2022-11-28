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


class Test_append_xr_to_nc:
    @pytest.fixture()
    def nc_dset(self):
        fname = uuid.uuid4()
        with nc.Dataset(filename=fname, mode='w', diskless=True) as dset:
            dset.createDimension('x', 4)
            dset.createDimension('y', None)
            dset.createVariable('a', 'i4', ('y', 'x'))
            dset.createVariable('b', 'i4', 'x')
            dset.variables['a'][:2, :] = 0
            dset.variables['b'][:] = 1
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


class Test_OutputNC:
    @pytest.fixture()
    def result(self):
        return xr.Dataset(
            data_vars=dict(
                x=xr.Variable('t', [0, 0]),
                y=xr.Variable('t', [0, 0]),
                z=xr.Variable('t', [0, 0]),
                u=xr.Variable('t', [0, 0]),
                v=xr.Variable('t', [0, 0]),
                w=xr.Variable('t', [0, 0]),
                density=xr.Variable('t', [0, 0]),
                radius=xr.Variable('t', [0, 0]),
            ),
            coords=dict(
                t=xr.Variable('t', [1000, 2000]),
            ),
        )

    def test_can_add_attributes(self, result):
        from uuid import uuid4
        with effluent.io.OutputNC(uuid4(), diskless=True) as out:
            out.write(time=0, result=result)
            assert out.dset.variables['z'].positive == 'down'

    def test_can_append_variables(self, result):
        from uuid import uuid4
        with effluent.io.OutputNC(uuid4(), diskless=True) as out:
            out.write(time=0, result=result)
            out.write(time=1, result=result)
            rt = out.dset.variables['release_time'][:]
            assert rt.tolist() == [0, 1]


class Test_OutputCSV:
    @pytest.fixture()
    def result(self):
        return xr.Dataset(
            data_vars=dict(
                x=xr.Variable('t', [0, 0]),
                y=xr.Variable('t', [0, 0]),
                z=xr.Variable('t', [0, 0]),
                u=xr.Variable('t', [0, 0]),
                v=xr.Variable('t', [0, 0]),
                w=xr.Variable('t', [0, 0]),
                density=xr.Variable('t', [0, 0]),
                radius=xr.Variable('t', [0, 0]),
            ),
            coords=dict(
                t=xr.Variable('t', [1000, 2000]),
            ),
        )

    def test_can_append_variables(self, result):
        from uuid import uuid4
        with effluent.io.OutputCSV(uuid4(), diskless=True) as out:
            out.write(time=0, result=result)
            out.write(time=1, result=result)

            out.dset.seek(0)
            lines = out.dset.readlines()

        assert len(lines) == 5
        assert lines[0] == 'release_time,t,x,y,z,u,v,w,density,radius\n'
        assert lines[1] == '0,1000,0,0,0,0,0,0,0,0\n'
        assert lines[2] == '0,2000,0,0,0,0,0,0,0,0\n'
        assert lines[3] == '1,1000,0,0,0,0,0,0,0,0\n'
        assert lines[4] == '1,2000,0,0,0,0,0,0,0,0\n'


class Test_Pipe_from_config:
    def test_explicit_mapping_with_arrays(self):
        conf = dict(
            time=[0, 1200], flow=[2, 3], dens=[4, 5], diam=[6, 7], depth=[8, 9],
            decline=[10, 11],
        )
        p = effluent.io.Pipe.from_config(conf)
        dset = p.select(time=600)
        assert dset.dens.values == 4.5

    def test_explicit_mapping_with_singlenums(self):
        conf = dict(
            time=[0, 1200], flow=[2, 3], dens=4, diam=6, depth=[8, 9],
            decline=10,
        )
        p = effluent.io.Pipe.from_config(conf)
        dset = p.select(time=600)
        assert dset.dens.values == 4

    def test_csv_file(self):
        buf = io.StringIO("""
            time, flow, dens, diam, depth, decline
               0,    1,    3,    5,     7,       9
            1200,    2,    4,    6,     8,      10
        """)
        conf = dict(csv=dict(file=buf))
        p = effluent.io.Pipe.from_config(conf)
        dset = p.select(time=600)
        assert dset.dens.values == 3.5

    def test_nc_file(self):
        buf = xr.Dataset(
            coords=dict(time=[0, 1200]),
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
        dset = p.select(time=600)
        assert dset.dens.values == 3.5


class Test_Ambient_from_config:
    def test_explicit_mapping_with_arrays(self):
        conf = dict(
            time=[0, 1200],
            depth=[0, 10, 20],
            coflow=[[0, 1, 2], [3, 4, 5]],
            crossflow=[[6, 7, 8], [9, 0, 1]],
            dens=[[2, 3, 4], [5, 6, 7]],
        )
        p = effluent.io.Ambient.from_config(conf)
        dset = p.select(time=600)
        assert dset.dens.values.tolist() == [3.5, 4.5, 5.5]

    def test_explicit_mapping_with_singlenums(self):
        conf = dict(
            time=[0, 1200],
            depth=[0, 10, 20],
            coflow=[[0, 1, 2], [3, 4, 5]],
            crossflow=[[6, 7, 8], [9, 0, 1]],
            dens=3,
        )
        p = effluent.io.Ambient.from_config(conf)
        dset = p.select(time=600)
        assert dset.dens.values.tolist() == [3, 3, 3]

    def test_csv_file(self):
        buf = io.StringIO("""
            time, depth, coflow, crossflow, dens
               0,     0,      0,         6,    2
               0,    10,      1,         7,    3
               0,    20,      2,         8,    4
            1200,     0,      3,         9,    5
            1200,    10,      4,         0,    6
            1200,    20,      5,         1,    7
        """)
        conf = dict(csv=dict(file=buf))
        p = effluent.io.Ambient.from_config(conf)
        dset = p.select(time=600)
        assert dset.dens.values.tolist() == [3.5, 4.5, 5.5]

    def test_nc_file(self):
        buf = xr.Dataset(
            coords=dict(time=[0, 1200], depth=[0, 10, 20]),
            data_vars=dict(
                coflow=xr.Variable(('time', 'depth'), [[0, 1, 2], [3, 4, 5]]),
                crossflow=xr.Variable(('time', 'depth'), [[6, 7, 8], [9, 0, 1]]),
                dens=xr.Variable(('time', 'depth'), [[2, 3, 4], [5, 6, 7]]),
            )
        ).to_netcdf()
        conf = dict(nc=dict(file=buf))
        p = effluent.io.Ambient.from_config(conf)
        dset = p.select(time=600)
        assert dset.dens.values.tolist() == [3.5, 4.5, 5.5]
