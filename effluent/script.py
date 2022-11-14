def run(*argv):
    conf_fname = argv[0]

    # Load config file
    import yaml
    with open(conf_fname, encoding='utf-8') as fp:
        conf = yaml.safe_load(fp)

    # Save output file
    import xarray as xr
    result_fname = conf['output']['file']
    xr.Dataset().to_netcdf(result_fname)

    return result_fname
