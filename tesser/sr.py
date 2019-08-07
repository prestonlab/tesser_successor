import numpy as np
import numpy.linalg as la
import scipy as sp
import scipy.spatial.distance as dist
import matplotlib.pyplot as plt
import learn
import util
import network


def make_envstep (DATA):
    ''' makes the envstep necessary for learn.py'''
    envstep = []
    objects = util.get_objects (DATA)
    index = 0

    for index in range ( len(objects)):
        obj = int(objects[index])
        try:
            envstep.append (obj)
        except:
            pass
    return envstep

def make_structured_data(SUBJECT):
    data, keys = util.read_files(SUBJECT, 'structured', [1, 2], list(range(1, 7)))
    runs = []
    for i in keys:
        part_num, run_num = i[-11:-10], i[-5:-4]
        runs.append([data[i], part_num, run_num])
    return runs


def explore_runs (OPTION, SUBJECT, GAMMA, ALPHA):

    SR_matrices = []
    part_run = []
    '''This loop adds the address, part number and run number to the runs array, so that the object
         sequence in each run can be inputted to the learning agent.
    '''
    runs = make_structured_data(SUBJECT)

    # This OPTION allows the SR matrix to persist across all runs from Part 1 and Part 2
    #     without ever resetting.
    if OPTION == 'persist':
        M = np.zeros([21,21])
        for run in runs:
            part_num, run_num = run[1], run[2]
            envstep = make_envstep (run[0])
            M = np.array(learn.run_experiment(envstep, GAMMA, ALPHA, np.copy(M)))
            SR_matrices.append (M)
            part_run.append([part_num, run_num])
            
    if OPTION == 'repeat':
        M = np.zeros([21,21])
        for time in range (100):
            for run in runs:
                part_num, run_num = run[1], run[2]
                envstep = make_envstep (run[0])
                M = np.array (learn.run_experiment (envstep, GAMMA, ALPHA, np.copy(M)))
        return M
    
    if OPTION == 'once':
        M = np.zeros([21,21])
        for run in runs:
            part_num, run_num = run[1], run[2]
            envstep = make_envstep (run[0])
            M = np.array (learn.run_experiment (envstep, GAMMA, ALPHA, np.copy(M)))
        return M
    
    # This OPTION allows the SR matrix to persist in Part 1 and Part 2, but resets it between them.
    if OPTION == 'reset':
        M = np.zeros ([21,21])
        is_reset = False
        for run in runs:
            part_num, run_num = run[1], run[2]
            if not is_reset and part_num == 2:
                M = np.zeros ([21,21])
                is_reset = True
            envstep = make_envstep (run[0])
            M = np.array(learn.run_experiment(envstep, GAMMA, ALPHA, np.copy(M)))
            SR_matrices.append (M)
            part_run.append([part_num, run_num])
    
    # This OPTION resets the SR matrix between each run.
    if OPTION == 'independent':
        for run in runs:
            part_num, run_num = run[1], run[2]
            M = np.zeros([21,21])
            envstep = make_envstep (run[0])
            M = learn.run_experiment(envstep, GAMMA, ALPHA, M)
            SR_matrices.append (M)
            part_run.append([part_num, run_num])

    # This OPTION forces the SR matrix to persist across all runs, but instead of plotting the SR matrix
    #     after each run, it plots the changes made to it after learning each object sequence.
    if OPTION == 'track changes':
        M = np.zeros([21,21])
        for run in runs:
            part_num, run_num = run[1], run[2]
            envstep = make_envstep (run[0])
            M_new = np.copy(M)
            M_new = learn.run_experiment(envstep, GAMMA, ALPHA, M_new)
            SR_matrices.append (M)
            part_run.append([part_num, run_num])
            M = M_new
            
    return SR_matrices, part_run

# Modify the OPTION in the following call to the main function in order to visualize the desired learning
#     sequence.
# explore_runs('track changes')

def compute_limit_matrix (gamma, adjacency):
    num_states = 21
    identity = np.eye (num_states)
    return np.linalg.inv (identity - gamma*adjacency/6)

def rda (matrix):
    return dist.squareform (dist.pdist (matrix, 'correlation'))

def correlate_rows (matrix):
    return np.dot(matrix, matrix.T) / (la.norm (matrix)**2)

def correlate_columns (matrix):
    return np.dot(matrix.T, matrix) / (la.norm (matrix)**2)

def pBGivenA (A, B, C, SR):
    return SR[A][B] / (SR[A][B] + SR[A][C])

def pCGivenA (A, B, C, SR):
    return 1 - pBGivenA (A, B, C, SR)

def likelihood (cue, opt1, opt2, response, SR):
    if response == 0:
        return pBGivenA (cue, opt1, opt2, SR)
    if response == 1:
        return pCGivenA (cue, opt1, opt2, SR)

def compute_correlations(OPTION, SUBJECT, GAMMA, ALPHA):
    nodes = network.node_info()
    adjacency = network.adjacency(nodes)
    transition = adjacency / 6
    L = compute_limit_matrix (.5, adjacency)
    L_vector = L.flatten()
    M = explore_runs ('once',SUBJECT, GAMMA, ALPHA)
    M_vector = M.flatten()
    plt.matshow (M, vmin = 0, vmax = .5)
    plt.colorbar()
    
    if OPTION == 'norm':
        print ('Norm of L - M: ')
        print (la.norm(L_vector - M_vector, np.inf))
    
    if OPTION == 'correlation':
        print ('Correlation of L, M: ')
        print (np.dot (L_vector, M_vector) / (la.norm (L_vector) * la.norm (M_vector)))
        
    plt.matshow (L, vmin = 0, vmax = .5)
    plt.colorbar()
    
# This function gives the probability of obtaining the choices in the run, given
#     specific values for alpha, gamma.
    
# For obj_sequence_directory, enter the folder where the runs, which give the object sequence,
#     are stored. Use tesser scan 100B. For induction_directory, enter the full path of the 
#     induction data.

def get_log_likelihood(SUBJECT, GAMMA, ALPHA):
    SR = explore_runs ('once', SUBJECT, GAMMA, ALPHA)
    data, keys = util.read_files(SUBJECT, 'induction')
    cue_sequence, opt1_sequence, opt2_sequence, response_sequence = util.get_induction_data (data[keys[0]])
    num_trials = len (cue_sequence)
    log_likelihood = 0
    for trial_num in range (num_trials):
        try:
            cue_num, opt1_num, opt2_num, response_num = int(cue_sequence[trial_num]) - 1, int(opt1_sequence[trial_num]) - 1, int(opt2_sequence[trial_num]) - 1, int(response_sequence[trial_num]) - 1
            trial_probability = likelihood (cue_num, opt1_num, opt2_num, response_num, SR)
            log_likelihood += np.log (trial_probability)
        except ValueError:
            pass
    
    return log_likelihood
