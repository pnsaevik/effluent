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

        self.start = 0
        self.stop = 0
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
        pipe = self.pipe.select(time)
        ambient = self.ambient.select(time)
        return pipe, ambient

    def run(self):
        for _ in self.irun():
            pass

    def irun(self):
        times = np.arange(self.start, self.stop + self.step / 2, self.step)

        with self.output as output:
            for time in times:
                self.solver.data = self.data(time)
                result = self.solver.solve()
                output.write(time, result)
                yield result


def load_config(fname_or_dict):
    """Load and parse input config, and convert to internal config format"""

    c = load_file_or_dict(fname_or_dict)

    release = c['output'].pop('release')
    solver = c.pop('solver', {})
    trajectory = c['output'].pop('trajectory')
    model = c.pop('model', {})

    # Convert dates to posix times
    np_epoch = np.datetime64('1970-01-01')
    for key in ['start', 'stop']:
        if key in release:
            dt64 = np.datetime64(release[key])
            posix = (dt64 - np_epoch) / np.timedelta64(1, 's')
            release[key] = posix

    c['timestepper'] = release
    c['solver'] = {**solver, **trajectory, **model}

    return c


def load_file_or_dict(f):
    if isinstance(f, dict):
        return f
    else:
        with open(f, 'rb') as fp:
            return toml.load(fp)
