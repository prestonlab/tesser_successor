# module that plots models in fit, sr and rsa
import numpy as np
import matplotlib.pyplot as plt
import sr
import fit


def plot_explore_runs(PATH, SUBJECT, OPTION, GAMMA, ALPHA):
    """ Program which creates plots for learning models in explore runs:
        INPUT:

        PATH: string describing the path taken to access tesser data
        SUBJECT: Integeger input representing a particular subject
        OPTION: String descrbing particular models to run. 
        ('persist', 'reset', 'independent', 'track', 'changes')
        GAMMA & ALPHA: discount and learning rate parameters. From 0.0 to 1.0.
    """
    figsize = (70, 30)
    size = 50
    fig, axs = plt.subplots(2, 6, figsize=figsize, sharex="col", sharey="row")
    plt.suptitle(
        "Learning: " + OPTION + "  SUBJECT: " + str(SUBJECT) + " with "
        r"$\gamma$ : " + str(GAMMA) + " and "
        r"$\alpha$ : " + str(ALPHA),
        size=80,
    )
    M, pr = sr.explore_runs(PATH, SUBJECT, OPTION, GAMMA, ALPHA)
    for i in range(len(M)):
        part = int(pr[i][0])
        run = int(pr[i][1])
        axs[part - 1, run - 1].matshow(M[i], cmap="viridis", vmin=0, vmax=1)
        axs[part - 1, run - 1].set_title(
            "Part_%s Run_%s" % (pr[i][0], pr[i][1]), size=size
        )

    fig.delaxes(axs[0][5])
    return plt.show()


def plot_adjecncy_matrix():
    nodes = network.node_info()
    adjacency = network.adjacency(nodes)
    transition = adjacency / 6
    L = sr.compute_limit_matrix(0.5, adjacency)
    plt.matshow(L, vmin=0, vmax=0.5)
    plt.colorbar()


def input_user(PATH, SUBJECT, TYPE, MODEL, GAMMA, ALPHA):
    if TYPE == "structured":
        
        plot_explore_runs(MODEL, SUBJECT, GAMMA, ALPHA, PATH)

    if TYPE == "induction":
        m1, m2 = fit.maximize_likelihood(SUBJECT)
        logl = fit.get_log_likelihood(SUBJECT, m2, m1)
        print(
            "The log likelihood for SUBJECT: %s is %s and is maximized with gamma: %s and alpha: %s"
            % (SUBJECT, logl, m2, m1)
        )

    else:
        print('Value Error TYPE == "structured" or "induction"')
