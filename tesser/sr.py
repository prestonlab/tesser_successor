# core model functions for learning/running experiment phases
import numpy as np
import numpy.linalg as la
from . import util
from . import network


class SR_Matrix:
    """ This class defines a reinforcement learning agent that
    learns the state-state successor representation without taking actions.
    Thus, the resulting SR matrix is in the service of prediction.
    Initialization parameters
    gamma: discount param
    alpha: learning rate
    p_sample: probability of sampling different options, only relevant for testing policy dependence
    NUM_STATES: the number of states in the environment to initialize matrices
    Ida Momennejad, 2019"""

    def __init__(self, gamma, alpha, NUM_STATES, M):
        self.gamma = gamma  # discount factor
        self.alpha = alpha  # learning rate
        self.M = M  # M: state-state SR
        self.onehot = np.eye(NUM_STATES)  # onehot matrix, for updating M

    def step(self, s, s_new):
        self.update_SR(s, s_new)

    def update_SR(self, s, s_new):
        self.M[s] = (1 - self.alpha) * self.M[s] + self.alpha * (
                self.onehot[s] + self.gamma * self.M[s_new]
        )


def run_experiment(envstep, gamma, alpha, M, N_STATES):
    """ This function uses the reinforcement learning agent class in
        SR_no_action.py to learn.
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
        NUM_STATES: the number of states in the environment to initialize matrices

        Outputs:
        M: SR matrix
        
        Ida Momennejad, NYC, 2019"""

    num_states = N_STATES

    SR_agent = SR_Matrix(gamma, alpha, num_states, M)
    s = envstep[0] - 1

    for entry in envstep[1:]:  # go through trajectory till the end

        s_new = entry - 1
        SR_agent.step(s, s_new)

        s = int(s_new)

    return SR_agent.M


def learn_sr(df, GAMMA, ALPHA):
    """Train an SR matrix on the structure-learning task."""

    n_states = len(np.unique(df.objnum))
    SR_matrices = {}
    M = np.zeros([n_states, n_states])
    for part in np.unique(df.part):
        for run in np.unique(df.loc[df.part == part, 'run']):
            envstep = df.loc[(df.part == part) & (df.run == run),
                             'objnum'].values
            M = np.array(run_experiment(envstep, GAMMA, ALPHA,
                                        np.copy(M), n_states))
            M = M / np.sum(M)
            SR_matrices[(part, run)] = M
    return SR_matrices


def explore_runs(df, OPTION, GAMMA, ALPHA):
    """This loop adds the address, part number and run number to the runs array, so that the object
         sequence in each run can be inputted to the learning agent.
        INPUT:

        df: Structured learning data in DataFrame form.
        OPTION: String describing particular models to run. 
        ('persist', 'repeat', 'once', 'reset', 'independent', 'track changes')
        GAMMA & ALPHA: discount and learning rate parameters. From 0.0 to 1.0.
    """

    n_states = len(np.unique(df.objnum))
    SR_matrices = {}
    data = []
    for part in (1, 2):
        runs = np.unique(df.loc[df.part == part, 'run'])
        for run in runs:
            # get data for this run
            df_run = df.loc[(df.part == part) & (df.run == run), :]
            obj = df_run.objnum.values
            data.append([obj, part, run])
    # This OPTION allows the SR matrix to persist across all runs from Part 1 and Part 2
    #     without ever resetting.
    if OPTION == "persist":
        M = np.zeros([n_states, n_states])
        for run in (range(0, 11)):
            part_num, run_num = data[run][1], data[run][2]
            envstep = data[run][0]
            M = np.array(run_experiment(envstep, GAMMA, ALPHA, np.copy(M), n_states))
            SR_matrices[(part_num, run_num)] = M

    if OPTION == "repeat":
        M = np.zeros([n_states, n_states])
        for time in range(100):
            for run in runs:
                envstep = data[run][0]
                M = np.array(run_experiment(envstep, GAMMA, ALPHA, np.copy(M), n_states))
        return M

    if OPTION == "once":
        M = np.zeros([n_states, n_states])
        for run in (range(0, 11)):
            envstep = data[run][0]
            M = np.array(run_experiment(envstep, GAMMA, ALPHA, np.copy(M), n_states))
        return M

    # This OPTION allows the SR matrix to persist in Part 1 and Part 2, but resets it between them.
    if OPTION == "reset":
        M = np.zeros([n_states, n_states])
        is_reset = False
        for run in (range(0, 11)):
            part_num, run_num = data[run][1], data[run][2]
            if not is_reset and part_num == 2:
                M = np.zeros([n_states, n_states])
                is_reset = True
            envstep = data[run][0]
            M = np.array(run_experiment(envstep, GAMMA, ALPHA, np.copy(M), n_states))
            SR_matrices[(part_num, run_num)] = M

    # This OPTION resets the SR matrix between each run.
    if OPTION == "independent":
        for run in (range(0, 11)):
            part_num, run_num = data[run][1], data[run][2]
            M = np.zeros([n_states, n_states])
            envstep = data[run][0]
            M = np.array(run_experiment(envstep, GAMMA, ALPHA, M, n_states))
            SR_matrices[(part_num, run_num)] = M

    # This OPTION forces the SR matrix to persist across all runs, but instead of plotting the SR matrix
    #     after each run, it plots the changes made to it after learning each object sequence.
    if OPTION == "track changes":
        M = np.zeros([n_states, n_states])
        for run in (range(0, 11)):
            part_num, run_num = data[run][1], data[run][2]
            envstep = data[run][0]
            M_new = np.copy(M)
            M_new = np.array(run_experiment(envstep, GAMMA, ALPHA, M_new, n_states))
            SR_matrices[(part_num, run_num)] = M_new - M
            M = M_new

    return SR_matrices


# Modify the OPTION in the following call to the main function in order to visualize the desired learning
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


def compute_correlations(DF, OPTION, GAMMA, ALPHA):
    """ Computes the norm or correlation between the SR matrix and the limit matrix, 
        given a subject and values for gamma, alpha.
        The norm option computes the infinity norm between the SR Matrix and the limit matrix.
        The correlation option computes the correlation between the SR Matrix and the limit matrix.
        INPUT:

        DF: Structured learning data in DataFrame format.
        OPTION: String describing particular function ( 'norm' or 'correlation'))
        GAMMA & ALPHA: discount and learning rate parameters. From 0.0 to 1.0.
    """
    n_states = len(np.unique(DF.objnum))
    nodes = network.node_info()
    adjacency = network.adjacency(nodes)
    L = compute_limit_matrix(0.5, adjacency, n_states)
    L_vector = L.flatten()
    M = explore_runs(DF, "once", GAMMA, ALPHA)
    M_vector = M.flatten()

    if OPTION == "norm":
        print("Norm of L - M: ")
        print(la.norm(L_vector - M_vector, np.inf))

    if OPTION == "correlation":
        print("Correlation of L, M: ")
        print(np.dot(L_vector, M_vector) / (la.norm(L_vector) * la.norm(M_vector)))
