import numpy as np
import logging
import tomli as toml

logger = logging.getLogger(__name__)


class Model:
    def __init__(self):
        self.pipe = None
        self.ambient = None
        self.output = None
        self.solver = None

        self.start = np.datetime64('1970-01-01')
        self.stop = np.datetime64('1970-01-01')
        self.step = 86400

    @staticmethod
    def from_config(fname_or_dict):
        from effluent.io import Pipe, Ambient, Output
        from effluent.solver import Solver

        conf = load_config(fname_or_dict)

        m = Model()
        m.pipe = Pipe.from_config(conf['pipe'])
        m.ambient = Ambient.from_config(conf['ambient'])
        m.output = Output.from_config(conf['output'])
        m.solver = Solver.from_config(conf['solver'])

        for key in ['start', 'stop', 'step']:
            if key in conf['timestepper']:
                value = conf['timestepper'][key]
                setattr(m, key, value)

        return m

    def data(self, time):
        pipe = self.pipe.select(time).compute()
        ambient = self.ambient.select(time).compute()
        return pipe, ambient

    def run(self):
        for _ in self.irun():
            pass

    def irun(self):
        # Compute times
        one_sec = np.timedelta64(1000000, 'us')
        duration = (np.datetime64(self.stop) - np.datetime64(self.start)) / one_sec
        num_times = int(duration / self.step) + 1
        times = np.datetime64(self.start) + np.arange(num_times) * one_sec * self.step

        # Do timestepping
        try:
            for time in times:
                self.solver.data = self.data(time)
                result = self.solver.solve()
                self.output.write(time, result)
                yield result

        finally:
            self.output.close()
            self.ambient.close()


def load_config(fname_or_dict):
    """Load and parse input config, and convert to internal config format"""

    c = load_file_or_dict(fname_or_dict)

    release = c['output'].pop('release')
    solver = c.pop('solver', {})
    trajectory = c['output'].pop('trajectory')
    model = c.pop('model', {})

    c['timestepper'] = release
    c['solver'] = {**solver, **trajectory, **model}

    return c


def load_file_or_dict(f):
    if isinstance(f, dict):
        return f
    else:
        with open(f, 'rb') as fp:
            return toml.load(fp)
