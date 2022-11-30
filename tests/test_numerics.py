import numpy as np
from effluent import numerics


class Test_bilin_inv:
    def test_can_invert_single_cell(self):
        x, y = numerics.bilin_inv(
            f=np.array([3, 3, 5, 5]),
            g=np.array([30, 50, 50, 30]),
            F=np.array([[0, 10], [0, 10]]),
            G=np.array([[0, 0], [100, 100]]),
        )
        assert x.tolist() == [0.3, 0.5, 0.5, 0.3]
        assert y.tolist() == [0.3, 0.3, 0.5, 0.5]

    def test_can_invert_when_integer_coordinate(self):
        f = np.array([0, 1, 0, 2])
        g = np.array([0, 0, 10, 10])
        F = np.array([[0, 1, 2], [0, 1, 2]])
        G = np.array([[0, 0, 0], [10, 10, 10]])
        x, y = numerics.bilin_inv(f, g, F, G)
        i = x.round().astype(int)
        j = y.round().astype(int)

        # Assert i and j are practically equal to x and y
        assert np.abs(i - x).max() < 1e-7
        assert np.abs(j - y).max() < 1e-7

        # Check that F, G interpolated to i, j gives f, g
        assert F[i, j].tolist() == f.tolist()
        assert G[i, j].tolist() == g.tolist()
