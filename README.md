# EasyAxisOrientationEstimation

This folder contains example code for the estimation of the spatial orientation of immobilized magnetic nanoparticles with parallel aligned easy axes. The method is sketched in the notebook provided here and described in detail in the associated publication.

To be announced

## How to get started

### Non-interactive Pluto Notebooks

The quickest way to get started is to explore the non-interactive snapshot of the notebook by opening the `orientationEstimation.html` in a browser, which illustrates the complete evaluation pipeline of the project. Note that the outputs of the notebooks are provided within this repository. They are described in more detail in the [Output](Output) section.

### Using Pluto Notebooks

In order to explore the code interactively, one first has to download and install the [Julia programming language](https://julialang.org/) (version 1.6 or later) as well as the [Pluto](https://github.com/fonsp/Pluto.jl#installation) notebook package (version 0.15 or later) from the julia REPL (read-eval-print loop):
```julia
using Pkg           # use package manager
Pkg.add("Pluto")    # add Pluto notebook package
```
Next, one has to [start up](https://github.com/fonsp/Pluto.jl#usage) a Pluto notebook server, which will direct you to a Pluto start page in your browser, from which you can start the `orientationEstimation.jl` notebook by navigating to the appropriate location.

### Without Pluto Notebooks

Alternatively, one can clone this repository and navigate to the folder in the command line and run the script from the julia REPL:
```julia
using Pkg           # use package manager
Pkg.activate(".")   # activate project environment
Pkg.instantiate()   # install packages required to run the script

# run script
include("orientationEstimation.jl")
```

## Output

The notebook contains the whole evaluation pipeline related to the project. Note that running the notebook for the first time can take serval minutes, as Julia packages are installed and a large number of data sets are processed.

Running the notebook/script will result in the (re)creation of
* an exemplary reconstruction (`images/eval1.gif`),
* a set of figures illustrating the estimation process (`eval2.gif`, `weights.svg`, `weightsandfit.svg`),
* the `results.csv`, which contains the easy axis orientation, its estimate, estimation error and the position of the sample inside the field of view, and
* a histogram with the distribution of the easy axis alignment estimation error (`errordistribution.svg`).

## Open MPI Data

The measurement data associated to this project is about 230 MB large and will be downloaded and stored automatically, when the code is executed for the first time.

The complete data set is published under a [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/legalcode) license and can be found here:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5336271.svg)](https://doi.org/10.5281/zenodo.5336271)
