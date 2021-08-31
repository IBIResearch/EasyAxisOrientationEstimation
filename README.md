# EasyAxisOrientationEstimation

This folder contains example code for the estimation of the spatial orientation of immobilized magnetic nanoparticles with parallel aligned easy axes. The method is sketched in the notebook provided here and described in detail in the associated publication.

To be announced

## How to get started

### Using Pluto Notebooks

In order to use this code one first has to download and install the [Julia programming language](https://julialang.org/) (version 1.6 or later) as well as the [Pluto](https://github.com/fonsp/Pluto.jl#installation) notebook package (version 0.15 or later). Next, one has to [start up](https://github.com/fonsp/Pluto.jl#usage) a Pluto notebook server, which will direct you to a Pluto start page in your browser, from which you can start the `orientationEstimation.jl` notebook by navigating to the appropriate location.

### Without Pluto Notebooks

Alternatively, one can clone this repository and navigate to the folder in the command line and run the script from the julia REPL (read-eval-print loop):
```julia
using Pkg           # use package manager
Pkg.activate(".")   # activate project environment
Pkg.instantiate()   # install packages required to run the script

# run script
include("orientationEstimation.jl")
```

### Non-interactive Pluto Notebooks

A non-interactive snapshot of the notebook can be explored by opening the `orientationEstimation.html` in a browser. 

## Output

The notebook contains the whole evaluation pipeline related to the project. Note that running the notebook for the first time can take serval minutes, as Julia packages are installed and a large number of data sets are processed.

Running the notebook will result in the (re)creation of
* an exemplary reconstruction (`images/eval1.gif`),
* a set of figures illustrating the estimation process (`eval2.gif`, `weights.svg`, `weightsandfit.svg`),
* the `results.csv`, which contains the easy axis orientation, its estimate, estimation error and the position of the sample inside the fov and
* a histogram with the distribution of the easy axis alignment estimation error (`errordistribution.svg`).

## Open MPI Data

The measurement data associated to this project is about 230 MB large and will be downloaded and stored automatically, when the code is executed for the first time.

The complete data set is published under a [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/legalcode) license and can be found here:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5336271.svg)](https://doi.org/10.5281/zenodo.5336271)
