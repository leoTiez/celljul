import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import cv2
import imageio


class Nucleation:
    def __init__(self, x, y, m, g, sig, t):
        self.x = int(x)
        self.y = int(y)
        self.m = m
        self.g = g
        self.sig = sig
        self.radius = 0.
        self.create_time = t

    def get_coordinates(self):
        return np.asarray([self.x, self.y])

    def update_radius(self, t):
        # fix radius growth to 1 per time step
        self.radius += self.g

    def calc_circular_radius(self):
        vol = self.sig * self.radius ** (self.m - 1)
        return int(np.sqrt(vol))


def simulation(t, state, time_state, is_dissociated, nucleation_list, m, n, g, sig, sim_protein=False, min_repair_t=3.):
    # Growth
    for nucleation in nucleation_list:
        nucleation.update_radius(t)
        circ_radius = nucleation.calc_circular_radius()
        state = cv2.circle(state, nucleation.get_coordinates(), circ_radius, 1., -1)

    mask = np.logical_and(state == 1., time_state == -1)
    time_state[mask] = t

    if sim_protein:
        is_dissociated = np.logical_or(
            np.logical_and(
                t - time_state > min_repair_t + np.random.exponential(scale=1.),
                time_state != -1.
            ),
            is_dissociated
        )
        state[is_dissociated] = 0

    nucleation_candidates = np.random.binomial(1, p=n, size=state.shape) == 1
    new_nucleation = zip(*np.where(np.logical_and(nucleation_candidates, time_state == -1.)))
    for x, y in new_nucleation:
        nucleation_list.append(Nucleation(x, y, m, g, sig, t))

    return state, nucleation_list


def run_simulation(
        m,
        theta,
        n_p,
        g_p,
        sig_p,
        to_time=80,
        shape=(200, 200),
        do_sort=False,
        save_gif=False,
        sim_protein=False,
        min_repair_t=3.
):
    state = np.zeros(shape)
    time_state = -np.ones(shape)
    is_dissociated = np.zeros(shape, dtype='bool')
    nucleation_list = []
    state_list = []
    sim_rep_frac = []

    plt.ion()
    norm = Normalize(vmin=0, vmax=1)
    plt.figure(figsize=(8, 7))
    ax = plt.gca()
    for t in np.arange(to_time):
        state, nucleation_list = simulation(
            t,
            state,
            time_state,
            is_dissociated,
            nucleation_list,
            m,
            n_p,
            g_p,
            sig_p,
            sim_protein,
            min_repair_t
        )
        ax.clear()
        if do_sort:
            y_sort = np.argsort(state).reshape(-1)
            x_sort = np.repeat(np.arange(state.shape[0]), state.shape[1]).reshape(-1)
            ax.imshow(state[x_sort, y_sort].reshape(state.shape), cmap='Reds', norm=norm)
            state_list.append(state[x_sort, y_sort].reshape(state.shape))
        else:
            ax.imshow(state, cmap='Reds', norm=norm)
            state_list.append(state.astype('uint8').copy() * 255)
        ax.set_title('Time %s' % t)
        fig = plt.gcf()
        fig.canvas.draw()

        coloured = np.sum(state == 1.)
        sim_rep_frac.append((coloured / state.size) * theta)

        fig.canvas.flush_events()

    if save_gif:
        save_name = 'sorted' if do_sort else 'random'
        imageio.mimsave('figures/gif/%s_500.gif' % save_name, state_list, format='GIF-PIL', fps=10)

    plt.ioff()
    plt.close('all')
    return sim_rep_frac

