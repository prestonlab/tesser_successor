# core model functions for learning/running experiment phases
import numpy as np
import numpy.linalg as la
from . import util
from . import network

# core model functions for learning/running experiment phases
import numpy as np


class SR_Matrix:
    """ This class defines a reinforcement learning agent that
    learns the state-state successor representation without taking actions.
    Thus, the resulting SR matrix is in the service of prediction.
    Initalization parameters
    gamma: discount param
    alpha: learning rate
    p_sample: probability of sampling different options, only relevant for testing poilcy dependence
    NUM_STATES: the number of states in the environment to intialize matrices
    Ida Momennejad, 2019"""

    def __init__(self, gamma, alpha, NUM_STATES, M):
        self.gamma = gamma  # discount factor
        self.alpha = alpha  # learning rate
        self.M = M  # M: state-state SR
        self.onehot = np.eye(NUM_STATES)  # onehot matrix, for updating M

    def step(self, s, s_new):
        self.update_SR(s, s_new)

    def update_SR(self, s, s_new):
        self.M[s] = (1 - self.alpha) * self.M[s] + self.alpha * (self.onehot[s] + self.gamma * self.M[s_new])


def run_experiment(envstep, gamma, alpha, M):
    """ This function uses the reinfrocement learning agent class in
        SR_no_action.py to learn.
        Here the function takes the environment from Experiment 1 in our
        Nat Hum Beh paper & learns predictive representations with the
        specified learning rate and scale.
        Note: This is not SR dyna, nor SR-MB.
        This agent only learns the SR.
        Inputs:

        envstep:  generated with ida_envs.generate_nathumbeh_env1()
        gamma: discount parameter, determines scale of predictive representations
        alpha: learning rate
        p_sample: prob sampling each of the two sequences
        verbose: TRUE or FALSE
        Outputs:
        M: SR matrix
        W: value weights W
        memory: memory of episodes
        episodies: # episodes it takes to reach convergence
        Ida Momennejad, NYC, 2019"""

    num_states = 21

    SR_agent = SR_Matrix(gamma, alpha, num_states, M)
    s = envstep[0] - 1

    for entry in envstep[1:]:  # go through trajectory till the end

        s_new = entry - 1
        SR_agent.step(s, s_new)

        s = int(s_new)

    return SR_agent.M




def make_envstep(DATA):
    """ makes the envstep necessary for 'learn' module"""
    envstep = []
    objects = util.get_objects(DATA)
    index = 0

    for index in range(len(objects)):
        obj = int(objects[index])
        try:
            envstep.append(obj)
        except:
            pass
    return envstep


def make_structured_data(PATH, SUBJECT):
    """ Returns a 3 column list to execute explore_runs models
        based on the structured learning data from tesser
        runs = [ Data for entry n, Part number for n, Run number for n ]
    """
    data, keys = util.read_files(PATH, SUBJECT, "structured", [1, 2], list(range(1, 7)))
    runs = []
    for i in keys:
        part_num, run_num = i[-11:-10], i[-5:-4]
        runs.append([data[i], part_num, run_num])
    return runs


def explore_runs(DF, OPTION, GAMMA, ALPHA):
    SR_matrices = []
    part_run = []
    num_runs = DF.shape[0]
    """This loop adds the address, part number and run number to the runs array, so that the object
         sequence in each run can be inputted to the learning agent.
        INPUT:

        PATH: string describing the path taken to access tesser data
        OPTION: String describing particular models to run. 
        ('persist', 'repeat', 'once', 'reset', 'independent', 'track', 'changes')
        SUBJECT: Integer input representing a particular subject
        GAMMA & ALPHA: discount and learning rate parameters. From 0.0 to 1.0.
    """
    runs = []
    for run_num in range (num_runs):
        run = []
        run.append (DF.at[run_num, 'Obj'])
        run.append (DF.at[run_num, 'PartNum'])
        run.append (DF.at[run_num, 'ObjNum'])
        runs.append (run)

    # This OPTION allows the SR matrix to persist across all runs from Part 1 and Part 2
    #     without ever resetting.
    if OPTION == "persist":
        M = np.zeros([21, 21])
        for run in runs:
            part_num, run_num = run[1], run[2]
            envstep = make_envstep(run[0])
            M = np.array(run_experiment(envstep, GAMMA, ALPHA, np.copy(M)))
            SR_matrices.append(M)
            part_run.append([part_num, run_num])

    if OPTION == "repeat":
        M = np.zeros([21, 21])
        for time in range(100):
            for run in runs:
                part_num, run_num = run[1], run[2]
                envstep = make_envstep(run[0])
                M = np.array(run_experiment(envstep, GAMMA, ALPHA, np.copy(M)))
        return M

    if OPTION == "once":
        M = np.zeros([21, 21])
        for run in runs:
            part_num, run_num = run[1], run[2]
            envstep = make_envstep(run[0])
            M = np.array(run_experiment(envstep, GAMMA, ALPHA, np.copy(M)))
        return M

    # This OPTION allows the SR matrix to persist in Part 1 and Part 2, but resets it between them.
    if OPTION == "reset":
        M = np.zeros([21, 21])
        is_reset = False
        for run in runs:
            part_num, run_num = run[1], run[2]
            if not is_reset and part_num == 2:
                M = np.zeros([21, 21])
                is_reset = True
            envstep = make_envstep(run[0])
            M = np.array(run_experiment(envstep, GAMMA, ALPHA, np.copy(M)))
            SR_matrices.append(M)
            part_run.append([part_num, run_num])

    # This OPTION resets the SR matrix between each run.
    if OPTION == "independent":
        for run in runs:
            part_num, run_num = run[1], run[2]
            M = np.zeros([21, 21])
            envstep = make_envstep(run[0])
            M = run_experiment(envstep, GAMMA, ALPHA, M)
            SR_matrices.append(M)
            part_run.append([part_num, run_num])

    # This OPTION forces the SR matrix to persist across all runs, but instead of plotting the SR matrix
    #     after each run, it plots the changes made to it after learning each object sequence.
    if OPTION == "track changes":
        M = np.zeros([21, 21])
        for run in runs:
            part_num, run_num = run[1], run[2]
            envstep = make_envstep(run[0])
            M_new = np.copy(M)
            M_new = run_experiment(envstep, GAMMA, ALPHA, M_new)
            SR_matrices.append(M)
            part_run.append([part_num, run_num])
            M = M_new

    return SR_matrices, part_run


# Modify the OPTION in the following call to the main function in order to visualize the desired learning
#     sequence.
# explore_runs('track changes')

def compute_limit_matrix(gamma, adjacency):
    """ Computes the matrix to which SR learning should converge, by summing a geometric matrix series."""
    num_states = 21
    identity = np.eye(num_states)
    return np.linalg.inv(identity - gamma * adjacency / 6)


def correlate_rows(matrix):
    """ Computes the correlation matrix for a matrix's rows."""
    return np.dot(matrix, matrix.T) / (la.norm(matrix) ** 2)


def correlate_columns(matrix):
    """ Computes the correlation matrix for a matrix's columns."""
    return np.dot(matrix.T, matrix) / (la.norm(matrix) ** 2)


def compute_correlations(PATH, SUBJECT, OPTION, GAMMA, ALPHA):
    """ Computes the norm or correlation between the SR matrix and the limit matrix, given a subject and values for gamma, alpha.
        The norm option computes the infinity norm between the SR Matrix and the limit matrix.
        The correlation option computes the correlation between the SR Matrix and the limit matrix.
        INPUT:

        PATH: string describing the path taken to access tesser data
        SUBJECT: Integer input representing a particular subject
        OPTION: String describing particular function ( 'norm' or 'correlation'))
        GAMMA & ALPHA: discount and learning rate parameters. From 0.0 to 1.0.
    """
    nodes = network.node_info()
    adjacency = network.adjacency(nodes)
    transition = adjacency / 6
    L = compute_limit_matrix(0.5, adjacency)
    L_vector = L.flatten()
    M = explore_runs(PATH, SUBJECT, "once", GAMMA, ALPHA)
    M_vector = M.flatten()

    if OPTION == "norm":
        print("Norm of L - M: ")
        print(la.norm(L_vector - M_vector, np.inf))

    if OPTION == "correlation":
        print("Correlation of L, M: ")
        print(np.dot(L_vector, M_vector) / (la.norm(L_vector) * la.norm(M_vector)))
