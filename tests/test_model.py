import xarray as xr
import numpy as np

from effluent import model
import pytest


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
