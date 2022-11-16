def run(*argv):
    conf_fname = argv[0]

    from .model import Model
    model = Model(conf_fname)
    model.run()

    return model.output.file
