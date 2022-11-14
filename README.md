# effluent

**NOTE**: This package is under development and not all features are
implemented yet.

Effluent is a python package for simulating the dispersion of effluent
discharges from wastewater pipes. The underlying model is based on
Lee, Joseph H. W., and Chu, Vincent H.: *Turbulent Jets and Plumes*.
Boston, MA: Springer US, 2003.
[doi:10.1007/978-1-4615-0407-8](https://doi.org/10.1007/978-1-4615-0407-8>).

The package is mainly intended for research purposes, and does not contain
any convenience plotting or statistics functionality. It is expected that
these analyses are conducted in post-processing stages using other packages.


## Installation

The package is installed using pip:

    pip install effluent
  

# Usage

The software is invoked from the command line as

    python -m effluent config.yaml

where `config.yaml` is the configuration file specifying the simulation
setup.


# Documentation

Check out the
[online documentation](https://effluent.readthedocs.io/en/latest/) for a
description of the algorithm and software features.
