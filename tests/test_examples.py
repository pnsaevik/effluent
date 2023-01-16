import pytest
import effluent
from pathlib import Path
import os
import tomli


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
        expected_fname = conf['output']['csv']['file']
        result_fname = 'result.csv'
        conf['output']['csv']['file'] = result_fname

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
    with open(a, 'r', encoding='utf-8') as fp:
        txt_a = fp.read()

    with open(b, 'r', encoding='utf-8') as fp:
        txt_b = fp.read()

    assert txt_a == txt_b


def compare_netcdf(a, b):
    import xarray as xr
    dset_a = xr.load_dataset(a)
    dset_b = xr.load_dataset(b)
    assert dset_a.to_dict() == dset_b.to_dict()


class Test_run:
    @pytest.mark.parametrize("name", NAMED_EXAMPLES)
    def test_matches_output(self, name):
        ex = Example(name)
        result, expected = ex.run()
        compare_files(result, expected)
        Path(result).unlink()
