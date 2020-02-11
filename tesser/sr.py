#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Learning module for tesser simulations. Functions for learning/running experiment phases
are based on Ida Momennejad's state-state successor representation learning agent.

These include:

- Class that defines a reinforcement learning agent which learns the state-state successor representation
without taking actions.
    SRMatrix()
- Function uses the reinforcement learning agent class in SRMatrix to learn.
    run_experiment(envstep, gamma, alpha, M, n_states)

- Train an SR matrix on the structure-learning task.
    learn_sr(df, gamma, alpha)

- Computes the matrix to which SR learning should converge, by summing a geometric matrix series.
    compute_limit_matrix(gamma, adjacency, n_states)

- Computes the correlation matrix for a matrix's rows.
    correlate_rows(matrix)

- Computes the correlation matrix for a matrix's columns.
    correlate_columns(matrix)

- Computes the norm or correlation between the SR matrix and the limit matrix, given a subject & values for gamma, alpha
    compute_correlations(df, option, gamma, alpha)

"""
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from matplotlib import colors
from . import network


class SRMatrix:
    """ This class defines a reinforcement learning agent that
    learns the state-state successor representation without taking actions.
    Thus, the resulting SR matrix is in the service of prediction.
    Initialization parameters
    gamma: discount param
    alpha: learning rate
    p_sample: probability of sampling different options, only relevant for testing policy dependence
    num_states: the number of states in the environment to initialize matrices
    Ida Momennejad, 2019"""

    def __init__(self, gamma, alpha, num_states, M):
        self.gamma = gamma  # discount factor
        self.alpha = alpha  # learning rate
        self.M = M  # M: state-state SR
        self.onehot = np.eye(num_states)  # onehot matrix, for updating M

    def step(self, s, s_new):
        self.update_SR(s, s_new)

    def update_SR(self, s, s_new):
        self.M[s] = (1 - self.alpha) * self.M[s] + self.alpha * (
                self.onehot[s] + self.gamma * self.M[s_new]
        )


def run_experiment(envstep, gamma, alpha, M, n_states):
    """ This function uses the reinforcement learning agent class in
        SRMatrix to learn.
        Here the function takes the environment from Experiment 1 in our
        Nat Hum Beh paper & learns predictive representations with the
        specified learning rate and scale.
        Note: This is not SR dyna, nor SR-MB.
        This agent only learns the SR.
        Inputs:

        envstep:  objects in a single run
        gamma: discount parameter, determines scale of predictive representations
        alpha: learning rate
        M: prob sampling each of the two sequences
        n_states: the number of states in the environment to initialize matrices

        Outputs:
        M: SR matrix
        
        Ida Momennejad, NYC, 2019"""

    num_states = n_states

    SR_agent = SRMatrix(gamma, alpha, num_states, M)
    s = envstep[0] - 1

    for entry in envstep[1:]:  # go through trajectory till the end

        s_new = entry - 1
        SR_agent.step(s, s_new)

        s = int(s_new)

    return SR_agent.M


def learn_sr(df, gamma, alpha):
    """Train an SR matrix on the structure-learning task."""

    n_states = len(np.unique(df.objnum))
    SR_matrices = {}
    M = np.zeros([n_states, n_states])
    for part in np.unique(df.part):
        for run in np.unique(df.loc[df.part == part, 'run']):
            envstep = df.loc[(df.part == part) & (df.run == run),
                             'objnum'].values
            M = np.array(run_experiment(envstep, gamma, alpha,
                                        np.copy(M), n_states))
            M = M / np.sum(M)
            SR_matrices[(part, run)] = M
    return SR_matrices


def explore_runs(df, option, gamma, alpha):
    """Function that takes in structured data and returns a dictionary or array that contains learned data.
        User has the ability to set learning parameters and how much data the learning agent can use in any given
        sequence.
        INPUT:
        df: Structured learning data in DataFrame form.
        option: String describing particular models to run. 
        ('repeat', 'once', 'reset', 'independent', 'track changes')
        gamma & alpha: discount and learning rate parameters. From 0.0 to 1.0.
         OUTPUT:
        For 'repeat' and 'once' an array will be given.
        For 'reset', 'independent' and 'track changes', a dictionary containing a DataFrame for
        each run and part combination given in the initial DataFrame.
    """

    n_states = len(np.unique(df.objnum))
    SR_matrices = {}
    M = np.zeros([n_states, n_states])

    if option == "repeat":
        for time in range(100):
            for part in np.unique(df.part):
                for run in np.unique(df.loc[df.part == part, 'run']):
                    envstep = df.loc[(df.part == part) & (df.run == run),
                                     'objnum'].values
                    M = np.array(run_experiment(envstep, gamma, alpha, np.copy(M), n_states))
                    M = M / np.sum(M)
        return M

    if option == "once":
        for part in np.unique(df.part):
            for run in np.unique(df.loc[df.part == part, 'run']):
                envstep = df.loc[(df.part == part) & (df.run == run),
                                 'objnum'].values
                M = np.array(run_experiment(envstep, gamma, alpha, np.copy(M), n_states))
                M = M / np.sum(M)
        return M

    # This option allows the SR matrix to persist in Part 1 and Part 2, but resets it between them.
    if option == "reset":
        for part in np.unique(df.part):
            if part == 2:
                M = np.zeros([n_states, n_states])
                for run in np.unique(df.loc[df.part == part, 'run']):
                    envstep = df.loc[(df.part == part) & (df.run == run),
                                     'objnum'].values
                    M = np.array(run_experiment(envstep, gamma, alpha, np.copy(M), n_states))
                    M = M / np.sum(M)
                    SR_matrices[(part, run)] = M
            else:
                for run in np.unique(df.loc[df.part == part, 'run']):
                    envstep = df.loc[(df.part == part) & (df.run == run),
                                     'objnum'].values
                    M = np.array(run_experiment(envstep, gamma, alpha, np.copy(M), n_states))
                    M = M / np.sum(M)
                    SR_matrices[(part, run)] = M

    # This option resets the SR matrix between each run.
    if option == "independent":
        for part in np.unique(df.part):
            for run in np.unique(df.loc[df.part == part, 'run']):
                M = np.zeros([n_states, n_states])
                envstep = df.loc[(df.part == part) & (df.run == run),
                                 'objnum'].values
                M = np.array(run_experiment(envstep, gamma, alpha, np.copy(M), n_states))
                M = M / np.sum(M)
                SR_matrices[(part, run)] = M

    # This option forces the SR matrix to persist across all runs, but instead of plotting the SR matrix
    #     after each run, it plots the changes made to it after learning each object sequence.
    if option == "track changes":
        for part in np.unique(df.part):
            for run in np.unique(df.loc[df.part == part, 'run']):
                envstep = df.loc[(df.part == part) & (df.run == run),
                                 'objnum'].values
                M_new = np.copy(M)
                M_new = np.array(run_experiment(envstep, gamma, alpha, M_new, n_states))
            
                SR_matrices[(part, run)] = M_new - M
                M = M_new

    return SR_matrices


# Modify the option in the following call to the main function in order to visualize the desired learning
#     sequence.
# explore_runs('track changes')


def compute_limit_matrix(gamma, adjacency, n_states):
    """ Computes the matrix to which SR learning should converge, by summing a geometric matrix series."""
    num_states = n_states
    identity = np.eye(num_states)
    return np.linalg.inv(identity - gamma * adjacency / 6)


def correlate_rows(matrix):
    """ Computes the correlation matrix for a matrix's rows."""
    return np.dot(matrix, matrix.T) / (la.norm(matrix) ** 2)


def correlate_columns(matrix):
    """ Computes the correlation matrix for a matrix's columns."""
    return np.dot(matrix.T, matrix) / (la.norm(matrix) ** 2)


def compute_correlations(df, option, gamma, alpha):
    """ Computes the norm or correlation between the SR matrix and the limit matrix, 
        given a subject and values for gamma, alpha.
        The norm option computes the infinity norm between the SR Matrix and the limit matrix.
        The correlation option computes the correlation between the SR Matrix and the limit matrix.
        INPUT:

        df: Structured learning data in DataFrame format.
        option: String describing particular function ( 'norm' or 'correlation'))
        gamma & alpha: discount and learning rate parameters. From 0.0 to 1.0.
    """
    n_states = len(np.unique(df.objnum))
    nodes = network.temp_node_info()
    adjacency = network.adjacency_mat(nodes)
    L = compute_limit_matrix(0.5, adjacency, n_states)
    L_vector = L.flatten()
    M = explore_runs(df, "once", gamma, alpha)
    M_vector = M.flatten()

    if option == "norm":
        print("Norm of L - M: ")
        print(la.norm(L_vector - M_vector, np.inf))

    if option == "correlation":
        print("Correlation of L, M: ")
        print(np.dot(L_vector, M_vector) / (la.norm(L_vector) * la.norm(M_vector)))


def plot_sr(SR, subject, option="Standard", gamma=0.5, alpha=0.5):
    fig, ax = plt.subplots(2, 6, figsize=(14, 6))
    plt.suptitle(
        "Learning: " + option + "  subject: " + str(subject) + " with "
                                                               "gamma : " + str(gamma) + " and "
                                                                                         "alpha : " + str(alpha)
    )
    images = []
    for i, part in enumerate((1, 2)):
        for j, run in enumerate(range(1, 7)):
            if (part, run) not in SR:
                fig.delaxes(ax[i, j])
                continue
            images.append(ax[i, j].matshow(SR[(part, run)]))
            ax[i, j].set_title("Part_%s Run_%s \n" % (part, run))
            ax[i, j].label_outer()

    vmin = min(image.get_array().min() for image in images)
    vmax = max(image.get_array().max() for image in images)
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    for im in images:
        im.set_norm(norm)
    cbar = fig.colorbar(images[1], ax=ax.ravel().tolist(), shrink=0.95)

    cbar.set_ticks(np.arange(0, 0.5, 0.01))
    # cbar.set_ticklabels()

    plt.show()
