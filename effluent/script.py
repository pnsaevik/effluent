def run(*argv):
    conf_fname = argv[0]

    init_logger()

    from .model import Model
    model = Model(conf_fname)
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
