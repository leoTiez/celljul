import numpy as np
import pymc3 as pm
from py.kjma import repair_fraction


def mcmc_sample(obs, time_points, m, theta, n_samples=2000, n_chains=5, v_type='argmin_n'):
    g = 1.  # fix growth rate of radius to 1 unit per time step
    with pm.Model():
        n = pm.Bound(pm.Exponential, lower=0.0, upper=1.0)('n', lam=1., testval=.01)
        sig = pm.Bound(pm.Normal, lower=0)('sig', mu=np.pi, sigma=1, testval=np.pi)
        beta = (sig * n * g**(m-1))**(1./m)
        p = pm.Deterministic('p', repair_fraction(time_points, m, beta, theta))
        obs = pm.Bernoulli('obs', p, observed=obs)

        step = pm.Metropolis()
        cellmd_trace = pm.sample(n_samples, step=step, chains=n_chains, return_inferencedata=False)

    n_p, g_p, sig_p = None, 1., None
    if 'argmin' in v_type.lower():
        if v_type.lower() == 'argmin_n':
            amin = np.argmin(cellmd_trace['n'])
        elif v_type.lower() == 'argmin_sig':
            amin = np.argmin(cellmd_trace['sig'])
        else:
            raise ValueError('Unsuppoted v_type')
        n_p = cellmd_trace['n'][amin]
        sig_p = cellmd_trace['sig'][amin]
    elif v_type.lower() == 'mean':
        n_p = np.mean(cellmd_trace['n'])
        sig_p = np.mean(cellmd_trace['sig'])
    else:
        raise ValueError('Unsupportede v_type')

    return n_p, g_p, sig_p
