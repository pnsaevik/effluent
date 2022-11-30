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
    def test_can_convert_multiple_twodim_data(self):
        # Create test data
        from io import StringIO
        buf = StringIO(
            'outer,inner,first,second\n'
            '200,100,1,10\n200,101,2,11\n200,102,3,12\n'
            '201,100,4,13\n201,101,5,14\n201,102,6,15\n'
        )
        buf.seek(0)

        # Convert to xarray.DataArray
        index_cols = ['outer', 'inner']
        dset = utils.csv_to_xr(buf, index_cols)

        assert list(dset.data_vars) == ['first', 'second']
        assert list(dset.coords) == ['outer', 'inner']
        assert dset.outer.values.tolist() == [200, 201]
        assert dset.inner.values.tolist() == [100, 101, 102]
        assert dset.first.values.tolist() == [[1, 2, 3], [4, 5, 6]]
        assert dset.second.values.tolist() == [[10, 11, 12], [13, 14, 15]]
