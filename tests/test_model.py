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
    def test_explicit_mapping_with_arrays(self):
        conf = dict(
            time=[0, 1], flow=[2, 3], dens=[4, 5], diam=[6, 7], depth=[8, 9],
            decline=[10, 11],
        )
        p = model.Pipe.from_config(conf)
        dset = p.select(time=0)
        assert dset.dens.values == 4


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
    @pytest.fixture(scope='class')
    def pipe_dset_horz(self):
        return xr.Dataset(dict(depth=100, u=1, w=0, dens=1000, diam=0.5))

    @pytest.fixture(scope='class')
    def pipe_dset_light(self):
        return xr.Dataset(dict(depth=100, u=1, w=0, dens=990, diam=0.5))

    @pytest.fixture(scope='class')
    def pipe_dset_decl(self):
        return xr.Dataset(dict(depth=100, u=1, w=1, dens=1000, diam=0.5))

    @pytest.fixture(scope='class')
    def ambient_dset_still(self):
        return xr.Dataset(
            data_vars=dict(
                u=xr.Variable('depth', [0, 0]),
                v=xr.Variable('depth', [0, 0]),
                dens=xr.Variable('depth', [1000, 1000]),
            ),
            coords=dict(depth=[0, 1000]),
        )

    @pytest.fixture(scope='class')
    def ambient_dset_cross(self):
        return xr.Dataset(
            data_vars=dict(
                u=xr.Variable('depth', [0, 0]),
                v=xr.Variable('depth', [1, 1]),
                dens=xr.Variable('depth', [1000, 1000]),
            ),
            coords=dict(depth=[0, 1000]),
        )

    @pytest.fixture(scope='class')
    def ambient_dset_coflow(self):
        return xr.Dataset(
            data_vars=dict(
                u=xr.Variable('depth', [0.5, 0.5]),
                v=xr.Variable('depth', [0, 0]),
                dens=xr.Variable('depth', [1000, 1000]),
            ),
            coords=dict(depth=[0, 1000]),
        )

    @pytest.fixture(scope='class')
    def ambient_dset_stratified(self):
        return xr.Dataset(
            data_vars=dict(
                u=xr.Variable('depth', [0, 0, 0]),
                v=xr.Variable('depth', [0, 0, 0]),
                dens=xr.Variable('depth', [980, 1000, 1020]),
            ),
            coords=dict(depth=[0, 100, 200]),
        )

    @pytest.fixture(scope='class')
    def result_horz_still(self, pipe_dset_horz, ambient_dset_still):
        steps = np.array([0, 10, 20])
        ivp = model.InitialValueProblem(steps, pipe_dset_horz, ambient_dset_still)
        return ivp.solve()

    @pytest.fixture(scope='class')
    def result_decl_still(self, pipe_dset_decl, ambient_dset_still):
        steps = np.array([0, 10, 20])
        ivp = model.InitialValueProblem(steps, pipe_dset_decl, ambient_dset_still)
        return ivp.solve()

    @pytest.fixture(scope='class')
    def result_horz_cross(self, pipe_dset_horz, ambient_dset_cross):
        steps = np.array([0, 10, 20])
        ivp = model.InitialValueProblem(steps, pipe_dset_horz, ambient_dset_cross)
        return ivp.solve()

    @pytest.fixture(scope='class')
    def result_horz_coflow(self, pipe_dset_horz, ambient_dset_coflow):
        steps = np.array([0, 10, 20])
        ivp = model.InitialValueProblem(steps, pipe_dset_horz, ambient_dset_coflow)
        return ivp.solve()

    @pytest.fixture(scope='class')
    def result_light_still(self, pipe_dset_light, ambient_dset_still):
        steps = np.array([0, 10, 20])
        ivp = model.InitialValueProblem(steps, pipe_dset_light, ambient_dset_still)
        return ivp.solve()

    @pytest.fixture(scope='class')
    def result_light_stratified(self, pipe_dset_light, ambient_dset_stratified):
        steps = np.linspace(0, 200, 11)
        ivp = model.InitialValueProblem(steps, pipe_dset_light, ambient_dset_stratified)
        return ivp.solve()

    @staticmethod
    def sign_of_change(dset, v):
        sgn = np.unique(np.sign(np.diff(dset[v].values)))
        if len(sgn) > 1:
            return np.nan
        else:
            return sgn[0]

    def test_returns_xr_dataset_when_horz_still(self, result_horz_still):
        assert isinstance(result_horz_still, xr.Dataset)

    def test_variables_change_as_expected_when_horz_still(self, result_horz_still):
        r = result_horz_still
        assert self.sign_of_change(r, 'x') == 1
        assert self.sign_of_change(r, 'y') == 0
        assert self.sign_of_change(r, 'z') == 0
        assert self.sign_of_change(r, 'u') == -1
        assert self.sign_of_change(r, 'v') == 0
        assert self.sign_of_change(r, 'w') == 0
        assert self.sign_of_change(r, 'density') == 0
        assert self.sign_of_change(r, 'radius') == 1

    def test_variables_change_as_expected_when_light_still(self, result_light_still):
        r = result_light_still
        assert self.sign_of_change(r, 'x') == 1
        assert self.sign_of_change(r, 'y') == 0
        assert self.sign_of_change(r, 'z') == -1
        assert self.sign_of_change(r, 'u') == -1
        assert self.sign_of_change(r, 'v') == 0
        assert self.sign_of_change(r, 'w') == -1
        assert self.sign_of_change(r, 'density') == 1
        assert self.sign_of_change(r, 'radius') == 1

    def test_variables_change_as_expected_when_horz_cross(self, result_horz_cross):
        r = result_horz_cross
        assert self.sign_of_change(r, 'x') == 1
        assert self.sign_of_change(r, 'y') == 1
        assert self.sign_of_change(r, 'z') == 0
        assert self.sign_of_change(r, 'u') == -1
        assert self.sign_of_change(r, 'v') == 1
        assert self.sign_of_change(r, 'w') == 0
        assert self.sign_of_change(r, 'density') == 0
        assert self.sign_of_change(r, 'radius') == 1

    def test_variables_change_as_expected_when_decl_still(self, result_decl_still):
        r = result_decl_still
        assert self.sign_of_change(r, 'x') == 1
        assert self.sign_of_change(r, 'y') == 0
        assert self.sign_of_change(r, 'z') == 1
        assert self.sign_of_change(r, 'u') == -1
        assert self.sign_of_change(r, 'v') == 0
        assert self.sign_of_change(r, 'w') == -1
        assert self.sign_of_change(r, 'density') == 0
        assert self.sign_of_change(r, 'radius') == 1

    def test_variables_change_as_expected_when_horz_coflow(self, result_horz_coflow):
        r = result_horz_coflow
        assert self.sign_of_change(r, 'x') == 1
        assert self.sign_of_change(r, 'y') == 0
        assert self.sign_of_change(r, 'z') == 0
        assert self.sign_of_change(r, 'u') == -1
        assert self.sign_of_change(r, 'v') == 0
        assert self.sign_of_change(r, 'w') == 0
        assert self.sign_of_change(r, 'density') == 0
        assert self.sign_of_change(r, 'radius') == 1

    def test_coflow_jet_travels_longer(self, result_horz_coflow, result_horz_still):
        x_coflow = result_horz_coflow.x.values
        x_still = result_horz_still.x.values
        assert x_coflow[-1] > x_still[-1]

    def test_variables_change_as_expected_when_light_stratified(self, result_light_stratified):
        r = result_light_stratified
        assert self.sign_of_change(r, 'x') == 1        # Increasing distance from outlet
        assert self.sign_of_change(r, 'y') == 0        # No cross-flow
        assert np.isnan(self.sign_of_change(r, 'z'))   # Vertical bouncing behaviour
        assert self.sign_of_change(r, 'u') == -1       # Decreasing horizontal speed
        assert self.sign_of_change(r, 'v') == 0        # No cross-flow
        assert np.isnan(self.sign_of_change(r, 'w'))   # Vertical bouncing behaviour
        assert np.isnan(self.sign_of_change(r, 'density'))  # Bouncing behaviour
        assert self.sign_of_change(r, 'radius') == 1   # Radius increase
