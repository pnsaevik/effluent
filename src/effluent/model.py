"""
The package contains the :class:`Model <effluent.model.Model>` class and helper
functions.
"""

import numpy as np
import logging
import tomli as toml
import effluent.io
import effluent.solver
import effluent.io
import effluent.solver
import xarray as xr

logger = logging.getLogger(__name__)


class Model:

    """
    A complete model of an effluent discharge.

    The class contains all parameters necessary to execute a simulation, including
    information about the pipe discharge itself, the ambient ocean, parameters for the
    mixing process, the desired output format and numerical solver parameters.

    :param pipe: Pipe and effluent discharge parameters
    :param ambient: Ambient ocean hydrography and currents
    :param output: Output target and output format
    :param solver: Physical parametrization and numerical solver
    :param start: First simulated discharge time
    :param stop: Last simulated discharge time
    :param step: Seconds between each simulated discharge time
    """

    def __init__(
            self, pipe: effluent.io.Pipe = None, ambient: effluent.io.Ambient = None,
            output: effluent.io.Output = None, solver: effluent.solver.Solver = None,
            start='1970-01-01', stop='1970-01-01', step=86400,
    ):
        self.pipe = pipe
        self.ambient = ambient
        self.output = output
        self.solver = solver

        self.start = start
        self.stop = stop
        self.step = step

    @staticmethod
    def from_config(fname_or_dict) -> "Model":
        """
        Initialize a model object using the :doc:`configuration format </config>`.

        :param fname_or_dict: Config dict or name of config file
        :return: An initialized object.
        """
        conf = load_config(fname_or_dict)

        m = Model()
        m.pipe = effluent.io.Pipe.from_config(conf['pipe'])
        m.ambient = effluent.io.Ambient.from_config(conf['ambient'])
        m.output = effluent.io.Output.from_config(conf['output'])
        m.solver = effluent.solver.Solver(**conf['solver'])

        for key in ['start', 'stop', 'step']:
            if key in conf['timestepper']:
                value = conf['timestepper'][key]
                setattr(m, key, value)

        return m

    def run(self):
        """
        Execute the simulation using repeated calls to
        :func:`self.irun() <effluent.model.Model.irun>`.

        The result is written to the output file.
        """
        for _ in self.irun():
            pass

    def irun(self) -> xr.Dataset:
        """
        Execute the simulation iteratively.

        For each iteration, a new discharge time is simulated using the following
        algorithm:

        1. Initial conditions are specified using
           :func:`self.solver.set_init() <effluent.solver.Solver.set_init>`

        2. Ambient conditions are specified using
           :func:`self.solver.set_ambient() <effluent.solver.Solver.set_ambient>`

        3. The numerical solution is computed using
           :func:`self.solver.solve() <effluent.solver.Solver.solve>`.

        4. The result is written to the output file using
           :func:`self.output.write() <effluent.io.Output.write>`

        5. The iterator yields the result

        :return: An iterator for simulating the discharge times sequentially.
        """

        # Compute times
        one_sec = np.timedelta64(1000000, 'us')
        duration = (np.datetime64(self.stop) - np.datetime64(self.start)) / one_sec
        num_times = int(duration / self.step) + 1
        times = np.datetime64(self.start) + np.arange(num_times) * one_sec * self.step

        # Do timestepping
        try:
            for time in times:
                self.solver.set_init(self.pipe, time)
                self.solver.set_ambient(self.ambient, time)
                result = self.solver.solve()
                self.output.write(time, result)
                yield result

        finally:
            self.output.close()
            self.ambient.close()


def load_config(fname_or_dict) -> dict:
    """
    Convert from external to internal config format.

    The :doc:`external config format </config>` is designed to be easy for the user,
    while the internal
    config format is structured as strictly one entry for each of the model
    subcomponent classes. The two formats are similar, but slightly reorganized.

    :param fname_or_dict: Either a dict, or the name of a file
    :return: A dict of config options, in internal config format
    """

    c = load_file_or_dict(fname_or_dict)

    release = c['output'].pop('release')
    solver = c.pop('solver', {})
    trajectory = c['output'].pop('trajectory')
    model = c.pop('model', {})

    c['timestepper'] = release
    c['solver'] = {**solver, **trajectory, **model}

    return c


def load_file_or_dict(f) -> dict:
    if isinstance(f, dict):
        return f
    else:
        with open(f, 'rb') as fp:
            return toml.load(fp)
