import pytest
import effluent
from pathlib import Path
import os


named_examples = ['minimal']


class Example:
    def __init__(self, name):
        self.name = name
        self.path = Path(__file__).parent / "examples" / name

    def run(self):
        os.chdir(self.path)
        conf_fname = self.path / 'config.yaml'
        result_fname = Path(effluent.run(conf_fname))
        expected_fname = self.path / 'expected.nc'
        if not expected_fname.exists():
            import shutil
            shutil.copyfile(result_fname, expected_fname)
        return result_fname, expected_fname


def compare_netcdf(a, b):
    import xarray as xr
    dset_a = xr.load_dataset(a)
    dset_b = xr.load_dataset(b)
    assert dset_a.to_dict() == dset_b.to_dict()


class Test_run:
    @pytest.mark.parametrize("name", named_examples)
    def test_matches_output(self, name):
        ex = Example(name)
        result, expected = ex.run()
        compare_netcdf(result, expected)
        result.unlink()
