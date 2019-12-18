# AnalySAS

[![Build Status](https://travis-ci.com/pozzo-reseach-group/AnalySAS.svg?branch=master)](https://travis-ci.com/pozzo-research-group/AnalySAS)

`analysas` is a Python package enabling common processing, analysis and modeling methods for small-angle scattering data. Some features include batch-removal of the incoherent background, efficient Porod analysis, and basic arithmetic transformations.

## Package Organization

The analysas package relies on the SasData class for handling scattering data. At the root level, one will find the basic functionalities of the package, e.g. arithmetic functions. Additional subpackages are described below.

### analysas.widgets

Includes pre-built interactive ipywidgets to visualize and/or perform data transformations, e.g. removing the incoherent background or trimming data points.

### analysas.models

Includes various models as well as methods to combine custom models for analyzing small-angle scattering data.

## Installation

Currently, the `analysas` package can be installed by downloading the contents of this repo, creating a Python environment with the package details found in either `environment.yml` or `requirements.txt`, and adding the path to AnalySAS to your system's search path when importing packages. Upcoming work will enable installation of the package via `pip`.

## Citing AnalySAS

We are working on finalizing proper citation of this work, but if you found the `analysas` package to be helpful in your own work, please reference the Pozzo Research Group at the University of Washington.
