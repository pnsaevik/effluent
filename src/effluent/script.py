def run(conf):
    """
    Run the main script and save the output in the specified file

    :param conf: Simulation configuration parameters (dict object or name of .toml file)
    See online documentation for a description of available options.

    :return: Name of output file
    """

    init_logger()

    from .model import Model
    model = Model.from_config(conf)
    model.run()

    return model.output.file


def init_logger(loglevel=None):
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
