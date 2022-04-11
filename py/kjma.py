import numpy as np


def repair_fraction(time, m, beta, theta):
    return (1 - np.exp(-(beta * time)**m)) * theta


def repair_fraction_over_time(to_time, m, beta, theta):
    return np.asarray([repair_fraction(t, m, beta, theta) for t in range(int(to_time))])


def create_data(m, beta, theta, to_time=120., n_points=100):
    obs_fraction = repair_fraction_over_time(to_time, m, beta, theta)
    # sample data points
    observation = np.asarray([
        np.random.binomial(1, p=np.minimum(np.maximum(p, 0), 1), size=n_points)
        for p in obs_fraction
    ]).reshape(-1)
    time_points = np.repeat(np.arange(to_time), n_points)
    return observation, time_points

