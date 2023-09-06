[![CircleCI](
https://circleci.com/gh/pnsaevik/effluent/tree/main.svg?style=shield)
](https://circleci.com/gh/pnsaevik/effluent)

# effluent

Effluent is a python package for simulating the dispersion of effluent
discharges from wastewater pipes. The underlying model is based on
Lee, Joseph H. W., and Chu, Vincent H.: *Turbulent Jets and Plumes*.
Boston, MA: Springer US, 2003.
[doi:10.1007/978-1-4615-0407-8](https://doi.org/10.1007/978-1-4615-0407-8>).

The package is mainly intended for research purposes, and does not contain
any convenience plotting or statistics functionality. It is expected that
these analyses are conducted in post-processing stages using other packages.


# Installation

The package can be installed using pip:

    pip install effluent

Alternatively, a [conda environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#) 
file (`environment.yml`) is included that creates an isolated python environment for Effluent.
The environment can be created by:

    conda env create -f environment.yml


# Usage

The software is invoked as a command line script:

    effluent config.toml

where `config.toml` is the configuration file specifying the simulation
setup. Examples of valid configuration files are given in the repository
directory `docs/source/examples`, and in the
[online documentation](https://effluent.readthedocs.io/en/latest/).


# Documentation

Check out the
[online documentation](https://effluent.readthedocs.io/en/latest/) for a
detailed description of the algorithm and software features.


# How to contribute

Use the GitHub issue tracking system to report problems with the software, seek
support, or suggest improvements. Code contributions can be suggested using
GitHub pull requests. Check documentation for details.
