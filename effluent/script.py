def run(*argv):
    conf_fname = argv[0]

    conf = load_config(conf_fname)
    pipe, ambient, output, solver = init_modules(conf)
    result = solver.solve(pipe, ambient)
    output.save(result)

    return output.file


def load_config(fname):
    import yaml
    with open(fname, encoding='utf-8') as fp:
        conf = yaml.safe_load(fp)

    if 'solver' not in conf:
        conf['solver'] = {}

    return conf


def init_modules(conf):
    pipe = Pipe(conf['pipe'])
    ambient = Ambient(conf['ambient'])
    output = Output(conf['output'])
    solver = Solver(conf['solver'])
    return pipe, ambient, output, solver


class Pipe:
    def __init__(self, conf):
        self.conf = conf


class Ambient:
    def __init__(self, conf):
        self.conf = conf


class Output:
    def __init__(self, conf):
        self.conf = conf
        self.file = conf['file']

    def save(self, result):
        from pathlib import Path
        suffix = Path(self.file).suffix

        if suffix == '.csv':
            with open(self.file, 'w', encoding='utf-8', newline='\n') as fp:
                from .utils import xr_to_csv
                xr_to_csv(result, fp)

        # Default file format: NetCDF4
        else:
            # darr.to_dataset().to_netcdf(fname)
            result.to_netcdf(self.file)


class Solver:
    def __init__(self, conf):
        self.conf = conf

    def solve(self, pipe, ambient):
        import xarray as xr
        return xr.Dataset()
