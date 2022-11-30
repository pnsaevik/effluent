from effluent import eos
import numpy as np
import pytest


@pytest.mark.parametrize("eos_name", ["roms"])
class Test_rho:
    @pytest.fixture()
    def rho(self, eos_name):
        return dict(roms=eos.roms_rho)[eos_name]

    def test_density_is_between_950_and_1100(self, rho):
        temp = np.array([0, 10, 20, 30, 40])
        salt = np.array([0, 40, 30, 0, 30])
        depth = np.array([0, 100, 10000, 500, 100])
        dens = rho(temp, salt, depth)
        assert np.all(dens > 950)
        assert np.all(dens < 1100)

    def test_density_increases_with_depth(self, rho):
        temp = 0
        salt = 0
        depth = np.linspace(0, 10000, 11)
        dens = rho(temp, salt, depth)
        assert np.all(np.diff(dens) > 0)

    def test_density_decreases_with_temp(self, rho):
        temp = np.linspace(4, 40, 11)
        salt = 0
        depth = 100
        dens = rho(temp, salt, depth)
        assert np.all(np.diff(dens) < 0)

    def test_density_increases_with_salt(self, rho):
        temp = 0
        salt = np.linspace(0, 50, 11)
        depth = 100
        dens = rho(temp, salt, depth)
        assert np.all(np.diff(dens) > 0)


class Test_roms_rho:
    def test_matches_check_value(self):
        dens = eos.roms_rho(temp=3, salt=35.5, depth=5000)
        expected = 1050.3639165364
        assert np.abs(dens - expected) < 1e-10
