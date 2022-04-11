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
I provide an interactive [Python Jupyer Notebook](KJMA\ Cell\ Dynamics.ipynb). Use `Python>=3`, I recommend `Python==3.6` or `Python==3.8`. Install requirements for the Jupyter Notebook with `pip`:

```commandline
python3 -m pip install -r requirements.txt
```

Initially, I intended to use this as an example to get a grip on Julia. However, I realised that the Python code seems to be more accurate. I provide therefore the code in both, Julia and Python. Use `Julia >= 1.`, I recommend `Julia==1.6`. To install the packages, run

```commandline
julia packages.jl
```
## Python Usage

I provide a command line module for Python. Run:

```commandline
python3 celljul.py -m=1.5 --beta=0.05 --theta=0.9 [-t=80 --dim=500 --sim_prot --do_sort --save_gif]
```

where parameters in brackets are optional. For more information type

```commandline
python3 celljul.py --help
```

## Julia Usage
Include the CellJul package via

```julia
include("celljul.jl)
```

and use the example as a reference: 
```julia
m = 1.5  # example m
beta = 1/20.  # example beta
theta = 0.9  # exmaple theta

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

n_p, g_p, sig_p = fetch_data(chain, v_type="argmin_n")

rep_frac = run_simulation(
    m, 
    beta,
    theta,
    n_p, 
    g_p,
    sig_p,
    sim_protein=true,  # Simulate the repair protein rather than the repair process
    to_time=120.,  # Simulation time
    dims=(100, 100),  # Dimensions of the grid. 
    verbosity=2,  # Verbosity level. Set it to 2 if you want to know the performance output
    to_gif=false  # Save simulation as a gif file
)
```