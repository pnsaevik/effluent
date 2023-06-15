"""
The module contains the main entry point of the application
"""


def run(conf):
    """
    Run the simulation and save the output in the specified file.

    Internally, the function creates a :class:`Model <effluent.model.Model>` object
    using the information from the config file, and calls the
    :func:`run <effluent.model.Model.run>` method of that object.

    :param conf: :doc:`Configuration parameters </config>` (dict object or file name).
    """

    init_logger()

    from .model import Model
    model = Model.from_config(conf)
    model.run()


def main():
    """
    Main script, runnable from the command line

    The function uses the `argparse` library to parse arguments from the command line.

    :return: 0 if successful
    """
    import argparse
    from . import __version__ as version_str

    parser = argparse.ArgumentParser(
        prog='effluent',
        description=(
            f'Effluent (v{version_str}) is a python package for simulating dispersion of effluent\n'
            'discharges from wastewater pipes.\n\n'
            'See online documentation at https://effluent.readthedocs.io/'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        'config_file',
        help='Path of config file',
    )

    args = parser.parse_args()

    run(args.config_file)

    return 0


def init_logger(loglevel=None):
    """
    Initialize the python logger

    :param loglevel: Desired output loglevel
    """

    import logging
    if loglevel is None:
        loglevel = logging.DEBUG

    package_name = str(__name__).split('.', maxsplit=1)[0]
    package_logger = logging.getLogger(package_name)
    package_logger.setLevel(loglevel)
    ch = logging.StreamHandler()
    ch.setLevel(loglevel)
    formatter = logging.Formatter('%(asctime)s  %(name)s:%(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    package_logger.addHandler(ch)
