# Effluent Tests

Running the tests requires adding `pytest` to the conda environment.
The file `test_environment.yml` includes this addition, so that `conda env create -f test_environment.yml` will create an appropriate environment for running the test suite.

To run the tests, simply activate the environment (`conda activate EffluentPyTest` and call `pytest`).
The tests should take around a minute to run.
