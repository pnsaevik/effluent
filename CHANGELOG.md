# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] - 2023-06-07

### Added
- Command line script for starting the program
- Citation guidance

### Fixed
- Included dask as an explicit dependency

## [1.2] - 2023-04-14
### Added
- Ambient input from ROMS ocean model
- Examples on how to use the software
- Adjustable float format in output CSV files
- Support for pandas 2.0
- In-repo paper for publication in Journal of Open Source Software

### Changed
- Use proper dates in input/output instead of seconds since simulation start

### Fixed
- Bug which caused command line invokation to fail unconditionally
- Correct computation of pipe flow speed from flow rate
- Computational elements are now oriented in the jet direction


## [1.1] - 2022-11-29
### Added
- Main script
- Configurable solver options
- Configurable model options
- Options for output to netCDF or CSV
- Options for input from netCDF, CSV or TOML
- Documented example
- Read The Docs documentation


## [1.0] - 2022-11-14

### Added

- Installable python package
