# RXN chemisttry utilities package

This repository contains various chemistry-related Python utilities used in the RXN universe.
For general utilities not related to chemistry, see our other repository [`rxn-utilities`](https://github.com/rxn4chemistry/rxn-utilities).

## TODO for open-sourcing

* GitHub actions
* Add CI badge
* Remove `.bumpversion.cfg`, `.whitesource`
* Rename to `rxn-chem-utils` in `setup.cfg`

## System Requirements

This package is supported on all operating systems. 
It has been tested on the following systems:
+ macOS: Big Sur (11.1)
+ Linux: Ubuntu 18.04.4

A Python version of 3.6 or greater is recommended.

## Installation guide

The package can be installed from Pypi:
```bash
pip install rxn-chem-utils
```

For local development, the package can be installed with:
```bash
pip install -e .[dev]
```

The `RDKit` dependency is not installed automatically and can be installed via Conda or Pypi:
```bash
# Install RDKit from Conda
conda install -c conda-forge rdkit

# Install RDKit from Pypi
pip install rdkit-pypi
```