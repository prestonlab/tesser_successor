# module that plots models in fit, sr and rsa
import numpy as np
import matplotlib.pyplot as plt
from . import sr
from . import fit
from . import network


def plot_explore_runs(SR, SUBJECT, OPTION, GAMMA, ALPHA):
    """ Program which creates plots for learning models in explore runs:
        INPUT:

        SR: Dictionary of SR matricies
        SUBJECT: Integeger input representing a particular subject
        OPTION: String descrbing particular models to run. 
        ('persist', 'reset', 'independent', 'track', 'changes')
        GAMMA & ALPHA: discount and learning rate parameters. From 0.0 to 1.0.
    """
    fig, ax = plt.subplots(2, 6, figsize=(14, 6))
    plt.suptitle(
        "Learning: " + OPTION + "  SUBJECT: " + str(SUBJECT) + " with "
        r"$\gamma$ : " + str(GAMMA) + " and "
        r"$\alpha$ : " + str(ALPHA)
    )
    for i, part in enumerate((1, 2)):
        for j, run in enumerate(range(1, 7)):
            if (part, run) not in SR:
                fig.delaxes(ax[i, j])
                continue
            im = ax[i,j].matshow(SR[(part, run)], vmin=0, vmax=.5)
            ax[i,j].set_title("Part_%s Run_%s \n" % (part, run))

    cbar = fig.colorbar(im, ax=ax.ravel().tolist(), shrink=0.95)

    cbar.set_ticks(np.arange(0, 0.5, 0.01))
    cbar.set_ticklabels(['low', 'medium', 'high'])

    plt.show()
    
def plot_rdms(rdms, SUBJECT,  GAMMA, ALPHA):
    """ Program which creates plots for learning models in explore runs:
        INPUT:


        SUBJECT: Integeger input representing a particular subject
        OPTION: String descrbing particular models to run. 
        ('persist', 'reset', 'independent', 'track', 'changes')
        GAMMA & ALPHA: discount and learning rate parameters. From 0.0 to 1.0.
    """
    fig, ax = plt.subplots(2, 6, figsize=(14, 6))
    plt.suptitle(
         "  SUBJECT: " + str(SUBJECT) + " with "
        r"$\gamma$ : " + str(GAMMA) + " and "
        r"$\alpha$ : " + str(ALPHA)
    )
    for i, part in enumerate((1, 2)):
        for j, run in enumerate(range(1, 7)):
            if (part, run) not in rdms:
                fig.delaxes(ax[i, j])
                continue
            im = ax[i,j].matshow(rdms[(part, run)], vmin=0, vmax=.5)
            ax[i,j].set_title("Part_%s Run_%s \n" % (part, run))

    cbar = fig.colorbar(im, ax=ax.ravel().tolist(), shrink=0.95)

    cbar.set_ticks(np.arange(0, 0.5, 0.01))
    cbar.set_ticklabels(['low', 'medium', 'high'])

    plt.show()

def plot_adjecncy_matrix():
    nodes = network.node_info()
    adjacency = network.adjacency(nodes)
    transition = adjacency / 6
    L = sr.compute_limit_matrix(0.5, adjacency, 21)
    plt.matshow(L, vmin=0, vmax=0.5)
    plt.colorbar()


# def input_user(DATAFRAME, SUBJECT, TYPE, MODEL, GAMMA, ALPHA):
#     if TYPE == "structured":
        
#         plot_explore_runs(MODEL, SUBJECT, GAMMA, ALPHA, PATH)

#     if TYPE == "induction":
#         m1, m2 = fit.maximize_likelihood(SUBJECT)
#         logl = fit.get_log_likelihood(SUBJECT, m2, m1)
#         print(
#             "The log likelihood for SUBJECT: %s is %s and is maximized with gamma: %s and alpha: %s"
#             % (SUBJECT, logl, m2, m1)
#         )

#     else:
#         print('Value Error TYPE == "structured" or "induction"')
