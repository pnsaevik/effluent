def run(*argv):
    import xarray as xr
    conf_fname = argv[0]

    conf = load_config(conf_fname)
    result = xr.Dataset()
    save_array(result, conf['output']['file'])

    return conf['output']['file']


def load_config(fname):
    import yaml
    with open(fname, encoding='utf-8') as fp:
        conf = yaml.safe_load(fp)

    return conf


def save_array(darr, fname):
    from pathlib import Path
    suffix = Path(fname).suffix

    if suffix == '.csv':
        with open(fname, 'w', encoding='utf-8', newline='\n') as fp:
            from .utils import xr_to_csv
            xr_to_csv(darr, fp)

    # Default file format: NetCDF4
    else:
        #darr.to_dataset().to_netcdf(fname)
        darr.to_netcdf(fname)
