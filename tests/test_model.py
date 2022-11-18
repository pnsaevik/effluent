import xarray as xr
import netCDF4 as nc
import numpy as np
from effluent import model
import pytest
import uuid


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
        model.write_xr_to_nc(xr_dset, nc_dset)
        assert nc_dset.variables['a'][:].tolist() == xr_dset['a'].values.tolist()

    def test_writes_data_coord_values(self, nc_dset):
        xr_dset = xr.Dataset(coords=dict(a=xr.Variable('a', np.arange(12))))
        model.write_xr_to_nc(xr_dset, nc_dset)
        assert nc_dset.variables['a'][:].tolist() == xr_dset['a'].values.tolist()

    def test_writes_variable_attrs(self, nc_dset):
        xr_dset = xr.Dataset(
            data_vars=dict(
                a=xr.Variable('x', np.arange(5), attrs=dict(units='m'))
            )
        )
        model.write_xr_to_nc(xr_dset, nc_dset)
        assert nc_dset.variables['a'].units == xr_dset['a'].attrs['units']

    def test_writes_dataset_attrs(self, nc_dset):
        xr_dset = xr.Dataset(
            data_vars=dict(a=xr.Variable('x', np.arange(5))),
            attrs=dict(date='2000-01-01'),
        )
        model.write_xr_to_nc(xr_dset, nc_dset)
        assert nc_dset.date == xr_dset.attrs['date']

    def test_writes_unlimited_dims(self, nc_dset):
        xr_dset = xr.Dataset(
            data_vars=dict(a=xr.Variable('x', np.arange(5))),
        )
        xr_dset.encoding['unlimited_dims'] = ['x']
        model.write_xr_to_nc(xr_dset, nc_dset)
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
        model.append_xr_to_nc(xr_dset, nc_dset)
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
        model.append_xr_to_nc(xr_dset, nc_dset)
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
        with model.OutputNC(uuid4(), diskless=True) as out:
            out.write(time=0, result=result)
            assert out.dset.variables['z'].positive == 'down'

    def test_can_append_variables(self, result):
        from uuid import uuid4
        with model.OutputNC(uuid4(), diskless=True) as out:
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
        with model.OutputCSV(uuid4(), diskless=True) as out:
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
    def test_can_interpret_mapping(self):
        conf = dict(
            format='dict',
            time=[0, 1],
            flow=[0, 0],
            dens=[0, 0],
            decline=0,
            diam=0,
            depth=0,
        )
        p = model.Pipe.from_config(conf)
        assert set(p.dset.variables) == {'u', 'w', 'diam', 'depth', 'dens', 'time'}


class Test_Ambient_from_config:
    def test_can_interpret_mapping(self):
        conf = dict(
            format='dict',
            time=[0, 1],
            depth=[0],
            coflow=[[0], [0]],
            crossflow=[[0], [0]],
            dens=[[0], [0]],
        )
        p = model.Ambient.from_config(conf)
        assert set(p.dset.variables) == {'u', 'v', 'depth', 'dens', 'time'}


class Test_IVP_solve:
    @pytest.fixture()
    def pipe_dset_horz(self):
        return xr.Dataset(dict(depth=100, u=1, w=0, dens=1000, diam=0.5))

    @pytest.fixture()
    def ambient_dset_still(self):
        return xr.Dataset(
            data_vars=dict(
                u=xr.Variable('depth', [0, 0]),
                v=xr.Variable('depth', [0, 0]),
                dens=xr.Variable('depth', [1000, 1000]),
            ),
            coords=dict(depth=[0, 1000]),
        )

    @pytest.fixture()
    def result_horz_still(self, pipe_dset_horz, ambient_dset_still):
        steps = np.array([0, 1000, 2000])
        ivp = model.InitialValueProblem(steps, pipe_dset_horz, ambient_dset_still)
        return ivp.solve()

    def test_returns_xr_dataset_when_horz_still(self, result_horz_still):
        assert isinstance(result_horz_still, xr.Dataset)

    def test_x_increases_when_horz_still(self, result_horz_still):
        x = result_horz_still.x.values
        assert np.all(np.diff(x) > 0)
