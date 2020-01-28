# module that plots models in fit, sr and rsa
import matplotlib.pyplot as plt
import numpy as np

from . import network
from . import sr


def plot_explore_runs(SR, subject, option, gamma, alpha, vmin=0, vmax=1):
    """ Program which creates plots for learning models in explore runs:
        INPUT:

        SR: Dictionary of SR matricies
        subject: Integeger input representing a particular subject
        option: String descrbing particular models to run. 
        ('persist', 'reset', 'independent', 'track changes')
        gamma & alpha: discount and learning rate parameters. From 0.0 to 1.0.
    """
    fig, ax = plt.subplots(2, 6, figsize=(14, 6))
    plt.suptitle(
        "Learning: " + option + "  subject: " + str(subject) + " with "
        r"$\gamma$ : " + str(gamma) + " and "
        r"$\alpha$ : " + str(alpha)
    )
    for i, part in enumerate((1, 2)):
        for j, run in enumerate(range(1, 7)):
            if (part, run) not in SR:
                fig.delaxes(ax[i, j])
                continue
            im = ax[i, j].matshow(SR[(part, run)], vmin=vmin, vmax=vmax)
            ax[i, j].set_title("Part_%s Run_%s \n" % (part, run))

    cbar = fig.colorbar(im, ax=ax.ravel().tolist(), shrink=0.95)

    cbar.set_ticks(np.arange(0, 0.5, 0.01))
    cbar.set_ticklabels(['low', 'medium', 'high'])

    plt.show()


def plot_rdms(rdms, subject, gamma, alpha):
    """ Program which creates plots for learning models in explore runs:
        INPUT:

        rdms: representational dissimilarity matrix for a list of matrices
        subject: Integeger input representing a particular subject
        gamma & alpha: discount and learning rate parameters. From 0.0 to 1.0.
    """
    fig, ax = plt.subplots(2, 6, figsize=(14, 6))
    plt.suptitle(
         "  subject: " + str(subject) + " with "
        r"$\gamma$ : " + str(gamma) + " and "
        r"$\alpha$ : " + str(alpha)
    )
    for i, part in enumerate((1, 2)):
        for j, run in enumerate(range(1, 7)):
            if (part, run) not in rdms:
                fig.delaxes(ax[i, j])
                continue
            im = ax[i, j].matshow(rdms[(part, run)], vmin=0, vmax=.5)
            ax[i, j].set_title("Part_%s Run_%s \n" % (part, run))

    cbar = fig.colorbar(im, ax=ax.ravel().tolist(), shrink=0.95)

    cbar.set_ticks(np.arange(0, 0.5, 0.01))
    cbar.set_ticklabels(['low', 'medium', 'high'])

    plt.show()


def plot_adjecncy_matrix():
    nodes = network.node_info()
    adjacency = network.adjacency(nodes)
    L = sr.compute_limit_matrix(0.5, adjacency, 21)
    plt.matshow(L, vmin=0, vmax=0.5)
    plt.colorbar()
