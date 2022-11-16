from effluent import utils
import numpy as np
import xarray as xr


class Test_xr_to_csv:
    def test_can_convert_arrays_of_twodim_data(self):
        # Create test dataset
        dset = xr.Dataset(
            data_vars=dict(
                first=xr.Variable(('outer', 'inner'), np.arange(1, 7).reshape((2, 3))),
                second=xr.Variable(('outer', 'inner'), np.arange(10, 16).reshape((2, 3))),
            ),
            coords=dict(inner=[100, 101, 102], outer=[200, 201]),
        )

        # Convert to csv
        from io import StringIO
        buf = StringIO()
        utils.xr_to_csv(dset, buf)

        # Test result
        assert buf.getvalue() == (
            'outer,inner,first,second\n'
            '200,100,1,10\n200,101,2,11\n200,102,3,12\n'
            '201,100,4,13\n201,101,5,14\n201,102,6,15\n'
        )


class Test_csv_to_xr:
    def test_can_convert_twodim_data(self):
        # Create test data
        from io import StringIO
        buf = StringIO(
            'outer,inner,mydata\n'
            '10,7,1\n10,8,2\n10,9,3\n'
            '11,7,4\n11,8,5\n11,9,6\n'
        )
        buf.seek(0)

        # Convert to xarray.DataArray
        darr = utils.csv_to_xr(buf)

        assert darr.name == 'mydata'
        assert darr.dims == ('outer', 'inner')
        assert darr.outer.values.tolist() == [10, 11]
        assert darr.inner.values.tolist() == [7, 8, 9]
        assert darr.values.tolist() == [[1, 2, 3], [4, 5, 6]]


class Test_bilin_inv:
    def test_can_invert_single_cell(self):
        x, y = utils.bilin_inv(
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
        x, y = utils.bilin_inv(f, g, F, G)
        i = x.round().astype(int)
        j = y.round().astype(int)

        # Assert i and j are practically equal to x and y
        assert np.abs(i - x).max() < 1e-7
        assert np.abs(j - y).max() < 1e-7

        # Check that F, G interpolated to i, j gives f, g
        assert F[i, j].tolist() == f.tolist()
        assert G[i, j].tolist() == g.tolist()
