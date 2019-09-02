# model fitting/parameter optimization
import numpy as np
from . import util
from . import sr
from scipy import optimize


def pBGivenA(A, B, C, SR):
    """ Computes the probability of B given A using the Luce choice rule."""
    return SR[A][B] / (SR[A][B] + SR[A][C])


def pCGivenA(A, B, C, SR):
    """ Computes the probability of C given A using the Luce choice rule."""
    return 1 - pBGivenA(A, B, C, SR)


def likelihood(cue, opt1, opt2, response, SR):
    """ Computes the likelihood of a subject's response given the SR matrix for this subject."""
    if response == 0:
        return pBGivenA(cue, opt1, opt2, SR)
    if response == 1:
        return pCGivenA(cue, opt1, opt2, SR)


def get_log_likelihood(PATH, SUBJECT, GAMMA, ALPHA):
    """ This function gives the probability of obtaining the choices in the run, given
        specific values for alpha, gamma.
        INPUT:

        PATH: string describing the path taken to access tesser data
        SUBJECT: Integeger input representing a particular subject
        GAMMA & ALPHA: discount and learning rate parameters. From 0.0 to 1.0.
    """
    SR = sr.explore_runs(PATH, SUBJECT, "once", GAMMA, ALPHA)
    data, keys = util.read_files(PATH, SUBJECT, "induction")
    cue_sequence, opt1_sequence, opt2_sequence, response_sequence = util.get_induction_data(
        data[keys[0]]
    )
    num_trials = len(cue_sequence)
    log_likelihood = 0
    for trial_num in range(num_trials):
        try:
            cue_num, opt1_num, opt2_num, response_num = (
                int(cue_sequence[trial_num]) - 1,
                int(opt1_sequence[trial_num]) - 1,
                int(opt2_sequence[trial_num]) - 1,
                int(response_sequence[trial_num]) - 1,
            )
            trial_probability = likelihood(
                cue_num, opt1_num, opt2_num, response_num, SR
            )
            log_likelihood += np.log(trial_probability)
        except ValueError:
            pass

    return log_likelihood


def maximize_likelihood(PATH, SUBJECT, OPTION):
    """ Numerically maximizes the log likelihood function on the set of the subject's choices
         to find optimal values for alpha, gamma
        INPUT:

        PATH: string describing the path taken to access tesser data
        SUBJECT: Integeger input representing a particular subject
    """
    def ll(x):
        ALPHA = x[0]
        GAMMA = x[1]
        return -get_log_likelihood (PATH, SUBJECT, ALPHA, GAMMA)
    
    if OPTION == 'basinhopping':
        alpha_max, gamma_max = optimize.basinhopping (ll, [1,1]).x
    if OPTION == 'differential evolution':
        alpha_max, gamma_max = optimize.differential_evolution (ll, [(0,1), (0,1)]).x
    if OPTION == 'brute':
        alpha_max, gamma_max = optimize.brute (ll, [(0,1), (0,1)])[0]
    
    return alpha_max, gamma_max
