#!/usr/bin/env python3
import sys
import argparse

import matplotlib.pyplot as plt

from py.kjma import *
from py.estimation import *
from py.simulation import *


def parse_cmd(args):
    parser = argparse.ArgumentParser(
        description='Visualise the KJMA transformation'
    )
    parser.add_argument('-m', type=float, required=True,
                        help='The m parameter of the KJMA model')
    parser.add_argument('--beta', type=float, required=True,
                        help='The beta parameter of the KJMA model')
    parser.add_argument('--theta', type=float, required=True,
                        help='The theta parameter of the KJMA model')
    parser.add_argument('-t', type=int, default=80,
                        help='Simulation time steps in minutes. Default is 80.')
    parser.add_argument('--sim_prot', action='store_true', dest='sim_prot',
                        help='If set, simulate association presence of abstract repair protein'
                             'rather than repair itself.')
    parser.add_argument('--do_sort', action='store_true', dest='do_sort',
                        help='If set, visualisation happens from right to left.')
    parser.add_argument('--save_gif', action='store_true', dest='save_gif',
                        help='If set, visualisation is saved as animation')
    parser.add_argument('--dim', type=int, default=500,
                        help='Quadratic size of simulated areas in positions')

    return parser.parse_args(args)


def main(p_args):
    m, beta, theta = p_args.m, p_args.beta, p_args.theta
    dim = p_args.dim
    to_time = p_args.t
    sim_prot = p_args.sim_prot
    do_sort = p_args.do_sort
    save_gif = p_args.save_gif
    obs, t_time = create_data(m, beta, theta, to_time=to_time)
    n_p, g_p, sig_p = mcmc_sample(obs, t_time, m, theta)
    sim_frac = run_simulation(
        m,
        theta,
        n_p,
        g_p,
        sig_p,
        shape=(dim, dim),
        sim_protein=sim_prot,
        do_sort=do_sort,
        save_gif=save_gif
    )

    plt.figure(figsize=(8, 7))
    plt.plot(repair_fraction_over_time(80, m, beta, theta), linestyle='--', linewidth=4, alpha=0.3, label='Model')
    plt.plot(sim_frac, linestyle='--', label='Simulation')
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel(r'$f(t)$', fontsize=21)
    plt.ylabel(r'Time (min)', fontsize=21)
    plt.legend(fontsize=21)
    plt.title('Model vs simulation', fontsize=32)
    plt.show()


if __name__ == '__main__':
    main(parse_cmd(sys.argv[1:]))


