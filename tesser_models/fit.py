# model fitting/parameter optimization
import numpy as np
import scipy.spatial.distance as spsd
from scipy.signal import find_peaks_cwt as find_peak
import numpy.linalg as la
import matplotlib.pyplot as plt
import util
import sr
import sys
sys.path.append ('/home/rodrigo/Documents/preston_labs/tesser_successor/')

from tesser import network

sr_data, dl, induc_data, idl = util.read_file()

def option (option, subj_num, alpha, gamma):
    figsize=(70,30)
    size = 50    
    fig, axs = plt.subplots(2, 6, figsize=figsize ,sharex='col', sharey='row')
    plt.subplots_adjust(hspace=0.5)
    a = util.initial_subj_run(subj_num)
    b = str(subj_num)

    # This option allows the SR matrix to persist across all runs from Part 1 and Part 2
    #     without ever resetting.
    if option == 'persist':
        plt.suptitle("Persistent Learning: Subject " + b + " Name: " + dl[a][:17] +  "\n", size = 80)
        M = np.zeros([21,21])
        for run in range(0,11): 
            if run < 5:
                M += sr.make_M_Matrix(a + run, gamma, alpha, np.copy(M))
                axs[0, run].matshow(M, cmap='viridis', vmin=0, vmax=1)
                axs[0, run].set_title(dl[a + run][-15:-4]+ "\n", size=size)
            else:    
                M += sr.make_M_Matrix(a + run, gamma, alpha, np.copy(M))
                axs[1, run -5].matshow(M, cmap='viridis', vmin=0, vmax=1)
                axs[1, run -5].set_title(dl[a + run][-15:-4]+ "\n", size=size)
        fig.delaxes(axs[0][5])
        return plt.show()

    # This option allows the SR matrix to persist in Part 1 and Part 2, but resets it between them.
    if option == 'reset':
        plt.suptitle("Reset Learning: Subject " + b + " Name: " + dl[a][:17] +  "\n", size = 80)
        M = np.zeros ([21,21])
        is_reset = False
        for run in range(0,11):
            if not is_reset and run >= 5:
                M = np.zeros ([21,21])
                is_reset = True
            if run < 5:
                M += sr.make_M_Matrix(a + run, gamma, alpha, np.copy(M))
                axs[0, run].matshow(M, cmap='viridis', vmin=0, vmax=1)
                axs[0, run].set_title(dl[a + run][-15:-4]+ "\n", size=size)
            else:    
                M += sr.make_M_Matrix(a + run, gamma, alpha, np.copy(M))
                axs[1, run -5].matshow(M, cmap='viridis', vmin=0, vmax=1)
                axs[1, run -5].set_title(dl[a + run][-15:-4]+ "\n", size=size)
        fig.delaxes(axs[0][5])
        return plt.show()

    # This option resets the SR matrix between each run.
    if option == 'independent':
        plt.suptitle("Independent Learning: Subject " + b + " Name: " + dl[a][:17] +  "\n", size = 80)
        M = np.zeros ([21,21])
        for run in range(0,11):
            if run < 5:
                M = sr.make_M_Matrix(a + run, gamma, alpha, M)
                axs[0, run].matshow(M, cmap='viridis', vmin=0, vmax=1)
                axs[0, run].set_title(dl[a + run][-15:-4]+ "\n", size=size)
                M = np.zeros ([21,21])
            else:    
                M = sr.make_M_Matrix(a + run, gamma, alpha, M)
                axs[1, run -5].matshow(M, cmap='viridis', vmin=0, vmax=1)
                axs[1, run -5].set_title(dl[a + run][-15:-4]+ "\n", size=size)
                M = np.zeros ([21,21])  
        fig.delaxes(axs[0][5])
        return plt.show()

    # This option forces the SR matrix to persist across all runs, but instead of plotting the SR matrix
    #     after each run, it plots the changes made to it after learning each object sequence.
    if option == 'track changes':
        plt.suptitle("Learning w/ Changes: Subject " + b + " Name: " + dl[a][:17] +  "\n", size = 80)
        M = np.zeros([21,21])
        for run in range(0,11):
            if run < 5:
                M_new = np.copy(M)
                M_new = sr.make_M_Matrix((a + run), gamma, alpha, np.copy(M))
                axs[0, run].matshow(M_new -M, cmap='viridis', vmin=0, vmax=1)
                axs[0, run].set_title(dl[a + run][-15:-4]+ "\n", size=size)
                M = M_new
            else:
                M_new = np.copy(M)
                M_new = sr.make_M_Matrix((a + run), gamma, alpha, np.copy(M))
                axs[1, run -5].matshow(M_new -M, cmap='viridis', vmin=0, vmax=1)
                axs[1, run -5].set_title(dl[a + run][-15:-4]+ "\n", size=size)
                M = M_new
        fig.delaxes(axs[0][5])
        return plt.show()

    if option == 'once':
        M = np.zeros([21,21])
        for run in range(0,11): 
            M += np.array ( sr.make_M_matrix((a + run), gamma, alpha, np.copy(M)))
        return M 

    if option == 'repeat':
        M = np.zeros([21,21])
        for time in range(0,101):
            for run in range(0,11):
                M = np.array ( sr.make_M_matrix((a + run), gamma, alpha, np.copy(M)))
        return M 


def indiv(subject_number, gamma, alpha):
    i = subject_number
    option('persist', i,gamma, alpha)
    option('reset', i,gamma, alpha)
    option('independent', i,gamma, alpha)
    option('track changes', i, gamma, alpha)

nodes = network.node_info()
adjacency = network.adjacency (nodes)
transition = adjacency / 6

def compute_limit_matrix (gamma):
    num_states = 21
    identity = np.eye (num_states)
    return np.linalg.inv (identity - gamma*adjacency/6)


def correlate_rows (matrix):
    return np.dot(matrix, matrix.T) / (la.norm (matrix)**2)

def correlate_columns (matrix):
    return np.dot(matrix.T, matrix) / (la.norm (matrix)**2)


def rain(option, subj_num, gamma, alpha):
    L = compute_limit_matrix (.5)
    L_vector = L.flatten()
    M = repeat(subj_num, gamma, alpha)
    M_vector = M.flatten()
#     plt.matshow (M, vmin = 0, vmax = .5)
#     plt.colorbar()
    
    if option == 'norm':
        opt = la.norm(L_vector - M_vector, np.inf)
        print ('Norm of L - M: ')
        print (opt)
    
    if option == 'correlation':
        opt = np.dot (L_vector, M_vector) / (la.norm (L_vector) * la.norm (M_vector))
#         print ('Correlation of L, M: ')
#         print (opt)

#     plt.matshow (L, vmin = 0, vmax = .5)
#     plt.colorbar()
    return opt


#from Demetriush Rowland
def pBGivenA (A, B, C, SR):
    return SR[A][B] / (SR[A][B] + SR[A][C])

def pCGivenA (A, B, C, SR):
    return 1 - pBGivenA (A, B, C, SR)

def likelihood (cue, opt1, opt2, response, SR):
    if response == 0:
        return pBGivenA (cue, opt1, opt2, SR)
    if response == 1:
        return pCGivenA (cue, opt1, opt2, SR)

def get_log_likelihood(subj_num, alpha, gamma):
    SR = option('once',subj_num, alpha, gamma)
    cue_sequence, opt1_sequence, opt2_sequence, response_sequence = sr.get_induc_data(subj_num)
    num_trials = len (cue_sequence)
    log_likelihood = 0
    for trial_num in range (num_trials):
        try:
            cue_num, opt1_num, opt2_num, response_num = int(cue_sequence[trial_num]) - 1, 
            int(opt1_sequence[trial_num]) - 1, int(opt2_sequence[trial_num]) - 1,
            int(response_sequence[trial_num]) - 1
            trial_probability = likelihood (cue_num, opt1_num, opt2_num, response_num, SR)
            log_likelihood += np.log (trial_probability)
        except ValueError:
            pass
    
    return log_likelihood


def main_log(subj_num):
    h = 10e-3
    alpha_max, gamma_max = 0.0, 0.0
    alpha, gamma = h, h
    log_likelihood_max = get_log_likelihood (subj_num, alpha, gamma)
    log_likelihoods = []
    
    while alpha <= 1:
        gamma = h
        while gamma <= 1:
            log_likelihood = get_log_likelihood(subj_num, alpha, gamma)
            log_likelihoods.append (log_likelihood)
            if log_likelihood > log_likelihood_max:
                alpha_max, gamma_max = alpha, gamma
            gamma += h
        alpha += h
    
    print(alpha_max, gamma_max)