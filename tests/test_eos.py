from effluent import eos
import numpy as np


class Test_rho:
    def test_density_is_between_950_and_1100(self):
        temp = np.array([0, 10, 20, 30, 40])
        salt = np.array([0, 40, 30, 0, 30])
        depth = np.array([0, 100, 10000, 500, 100])
        dens = eos.rho(temp, salt, depth)
        assert np.all(dens > 950)
        assert np.all(dens < 1100)

    def test_density_increases_with_depth(self):
        temp = 0
        salt = 0
        depth = np.linspace(0, 10000, 11)
        dens = eos.rho(temp, salt, depth)
        assert np.all(np.diff(dens) > 0)

    def test_density_decreases_with_temp(self):
        temp = np.linspace(0, 50, 11)
        salt = 0
        depth = 100
        dens = eos.rho(temp, salt, depth)
        assert np.all(np.diff(dens) < 0)

    def test_density_increases_with_salt(self):
        temp = 0
        salt = np.linspace(0, 50, 11)
        depth = 100
        dens = eos.rho(temp, salt, depth)
        assert np.all(np.diff(dens) > 0)
