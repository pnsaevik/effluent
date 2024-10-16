import pandas as pd
import pytest
import effluent
from pathlib import Path
import os
import tomli
import numpy as np

# noinspection PyUnresolvedReferences
import netCDF4 as nc  # Must be imported together with numpy to avoid warning


EXAMPLE_DIR = Path(__file__).parent.parent / "docs/source/examples"
NAMED_EXAMPLES = [f.name for f in EXAMPLE_DIR.glob('*') if f.is_dir()]


class Example:
    def __init__(self, name):
        self.name = name
        self.path = EXAMPLE_DIR / name

    def run(self):
        # Load and modify config file
        os.chdir(self.path)
        conf_fname = self.path / 'config.toml'
        with open(conf_fname, 'br') as fp:
            conf = tomli.load(fp)

        if 'csv' in conf['output']:
            output_type = 'csv'
        else:
            output_type = 'nc'

        expected_fname = conf['output'][output_type]['file']
        result_fname = 'result.' + output_type
        conf['output'][output_type]['file'] = result_fname

        # Run simulation
        effluent.run(conf)

        # For convenience: If "expected" file is nonexistent, create it
        if not Path(expected_fname).exists():
            import shutil
            shutil.copyfile(result_fname, expected_fname)

        return result_fname, expected_fname


def compare_files(a, b):
    suffix = Path(a).suffix
    comparators = {'.nc': compare_netcdf, '.csv': compare_csv}
    comparator = comparators[suffix]
    comparator(a, b)


def compare_csv(a, b):
    df_a = pd.read_csv(a)
    df_b = pd.read_csv(b)

    assert list(df_a.columns) == list(df_b.columns)
    for colname in df_a.columns:
        vals_a = df_a[colname].values
        vals_b = df_b[colname].values
        assert smalldiff(vals_a, vals_b, tol=1e-4)


def smalldiff(a, b, tol=1e-7):
    if np.issubdtype(a.dtype, np.integer):
        return a.tolist() == b.tolist()
    elif np.issubdtype(a.dtype, np.datetime64):
        return a.tolist() == b.tolist()
    elif isinstance(a[0], str):
        return list(a) == list(b)
    else:
        scale = np.max(tol + np.abs(a) + np.abs(b))
        diff = (a - b) / scale
        return np.max(np.abs(diff)) < tol


def compare_netcdf(a, b):
    import xarray as xr
    dset_a = xr.load_dataset(a)
    dset_b = xr.load_dataset(b)

    assert dset_a.attrs == dset_b.attrs
    assert set(dset_a.variables.keys()) == set(dset_b.variables.keys())
    for k in dset_a.variables.keys():
        assert dset_a.variables[k].dims == dset_b.variables[k].dims
        assert dset_a.variables[k].attrs == dset_b.variables[k].attrs
        assert dset_a.variables[k].shape == dset_b.variables[k].shape
        assert smalldiff(dset_a.variables[k].values, dset_b.variables[k].values)


class Test_run:
    @pytest.mark.parametrize("name", NAMED_EXAMPLES)
    def test_matches_output(self, name):
        ex = Example(name)
        result, expected = ex.run()
        compare_files(result, expected)
        Path(result).unlink()
