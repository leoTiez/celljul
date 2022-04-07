# CellJul

This repository is an extension of a previous study that modelled repair processes (or any cellular
process that can be understood as a system with a binary state transition) using the KJMA model (see
[GitHub repository](https://github.com/leoTiez/jmak) and publication on [BioRxiv](https://doi.org/10.1101/2022.03.29.486283)). It intends to simulate and visualise the processes that are predicted with the
described model. A two-dimensional understanding can provide the missing link to other dynamic processes,
such as the movement of repair proteins. The following two gifs show a simulation of the expansion of repair positions along a two-dimensional cell culture (left) as well as the movement of an abstract repair protein.

| | |
:---:|:---:
 ![](figures/gif/simulation_repair_(100,%20100).gif) | ![](figures/gif/simulation_protein_(100,%20100).gif) |

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

We recommend to use the default setting, which defines volume growth to be 1 per time step. If you know the KJMA parameters $m$, $\beta$, and $\theta$, you can simply run 

```julia
m = 1.5  # example m
beta = 1/20.  # example beta
theta = 0.9  # exmaple theta

rep_frac = run_simulation(
    m, 
    beta,
    theta,
    sim_protein=true,  # Simulate the repair protein rather than the repair process
    to_time=120.,  # Simulation time
    dims=(100, 100),  # Dimensions of the grid. 
    verbosity=2,  # Verbosity level. Set it to 2 if you want to know the performance output
    to_gif=false  # Save simulation as a gif file
)
```

If you want to determine the full range of possible values for the nucleation rate $n$, the growth speed $g$, and the shape $\sigma$, you can run the approximation method. This requires the creation of mock-up data.

```julia
obs, t_time = create_data(
    120.,  # time scale for the mock up data
    params,  # KJMA parameters as a vector
    n_points=100  # Number of data points per time point
)

chain = mcmc_sample(
    obs,
    t_time, 
    m, 
    theta,
    exp_decay=1.,  # Exponential decay for the sampling of n and g
    norm_sig=1.,  # Variance of the shape parameter sigma  
    iterations=1000,  # Number of iterations
    show_progress=true  # Show progress bar
)

n, g, m = fetch_data(chain, v_type="argmin_n")
```