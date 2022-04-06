# CellJul

This repository is an extension of a previous study that modelled repair processes (or any cellular
process that can be understood as a system with a binary state transition) using the KJMA model (see
[GitHub repository](https://github.com/leoTiez/jmak) and publication on [BioRxiv](https://doi.org/10.1101/2022.03.29.486283)). It intends to simulate and visualise the processes that are predicted with the
described model. A two-dimensional understanding can provide the missing link to other dynamic processes,
such as the movement of repair proteins. The following gif shows a simulation of the expansion of repair positions along a two-dimensional cell culture.

![CellJul Gif](figures/gif/random_500.gif)

## Install
We provide an interactive [Python Jupyer Notebook](KJMA\ Cell\ Dynamics.ipynb) (however, at the curren stage, there are no explanations). Use `Python>=3`, we recommend `Python==3.6` or `Python==3.8`. Install requirements for the Jupyter Notebook with `pip`:

```commandline
python3 -m pip install -r requirements.txt
```

The rest of the code is implemented in Julia. Use `Julia >= 1.`, we recommend `Julia==1.6`. To install the packages, run

```commandline
julia packages.jl
```

## Usage
Include the CellJul package via

```julia
include("celljul.jl)
```

We will shortly provide an example case.