import numpy as np
import logging

from effluent.solver import Solver

logger = logging.getLogger(__name__)


class Model:
    def __init__(self, fname_or_dict):
        from effluent.io import load_config, Pipe, Ambient, Output

        self.conf = load_config(fname_or_dict)
        self.pipe = Pipe.from_config(self.conf['pipe'])
        self.ambient = Ambient.from_config(self.conf['ambient'])
        self.output = Output.from_config(self.conf['output'])
        self.solver = Solver.from_config(self.conf['solver'])

    def data(self, time):
        pipe = self.pipe.select(time)
        ambient = self.ambient.select(time)
        return pipe, ambient

    def run(self):
        frequency = self.conf['timestepper']['frequency']
        stop = self.conf['timestepper']['stop']
        times = np.arange(0, stop + frequency / 2, frequency)

        with self.output as output:
            for time in times:
                self.solver.data = self.data(time)
                result = self.solver.solve()
                output.write(time, result)
