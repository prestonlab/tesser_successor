# model fitting/parameter optimization
import numpy as np
from . import util
from . import sr
from scipy.spatial import distance
from scipy import optimize
import time


def eu_dist(a, b, SR):
    return distance.euclidean(SR[a], SR[b])

# def pBGivenA(A, B, C, SR, tau):
#     """ Computes the probability of B given A using the Luce choice rule."""
#     if SR[A, B] == 0 and SR[A, C] == 0:
#         # if SR is zero for both, probability is undefined
#         return np.nan
#     return np.exp((-eu_dist(A,B,SR)/tau)) / (np.exp((-eu_dist(A,B,SR)/tau)) + np.exp((-eu_dist(A,C,SR)/tau)))

def pBGivenA(A, B, C, SR, tau):
    """ Computes the probability of B given A using the Luce choice rule."""
    if SR[A, B] == 0 and SR[A, C] == 0:
        # if SR is zero for both, probability is undefined
        return np.nan
    return SR[A][B]**tau / (SR[A][B]**tau + SR[A][C]**tau)


def pCGivenA(A, B, C, SR, tau):
    """ Computes the probability of C given A using the Luce choice rule."""
    return 1 - pBGivenA(A, B, C, SR, tau)


def likelihood(cue, opt1, opt2, response, SR, tau):
    """ Computes the likelihood of a subject's response given the SR matrix for this subject."""
    if response == 0:
        return pBGivenA(cue, opt1, opt2, SR, tau)
    if response == 1:
        return pCGivenA(cue, opt1, opt2, SR, tau)



def get_log_likelihood(STRUC_DF, INDUC_DF, GAMMA, ALPHA, tau, RETURN_TRIAL=False):
    """ This function gives the probability of obtaining the choices in the run, given
        specific values for alpha, gamma.
        INPUT:

        STRUC_DF: Structured learning data in DataFrame format.
        INDUC_DF: Generalized induction data in DataFrame format.
        GAMMA & ALPHA: discount and learning rate parameters. From 0.0 to 1.0.
    """
    SR = sr.explore_runs(STRUC_DF, "once", GAMMA, ALPHA)
    cue_sequence, opt1_sequence, opt2_sequence, response_sequence = util.get_induction_data(
        INDUC_DF
    )
    num_trials = len(cue_sequence)
    log_likelihood = 0
    all_trial_prob = np.zeros(num_trials)
    for trial_num in range(num_trials):
        if np.isnan(response_sequence[trial_num]):
            trial_probability = 0.5
            log_likelihood += np.log(trial_probability)
            all_trial_prob[trial_num] = trial_probability
            continue
        cue_num, opt1_num, opt2_num, response_num = (
            int(cue_sequence[trial_num]) - 1,
            int(opt1_sequence[trial_num]) - 1,
            int(opt2_sequence[trial_num]) - 1,
            int(response_sequence[trial_num]) - 1,
        )
        trial_probability = likelihood(
            cue_num, opt1_num, opt2_num, response_num, SR, tau
        )
        eps = 0.000001
        if np.isnan(trial_probability):
            # probability undefined; can occur if SR has zeros
            trial_probability = eps
        elif trial_probability < eps:
            # probability too close to zero; set to minimal value
            trial_probability = eps
        log_likelihood += np.log(trial_probability)
        all_trial_prob[trial_num] = trial_probability

    if RETURN_TRIAL:
        return log_likelihood, all_trial_prob
    else:
        return log_likelihood


def maximize_likelihood(STRUC_DF, INDUC_DF, OPTION,tau):
    """ Numerically maximizes the log likelihood function on the set of the subject's choices
         to find optimal values for alpha, gamma
        INPUT:

        STRUC_DF: Structured learning data in DataFrame format.
        INDUC_DF: Generalized induction data in DataFrame format.
        OPTION: scipy optimization functions:'basinhopping','differential evolution', 'brute'
    """
    def ll(x):
        ALPHA = x[0]
        GAMMA = x[1]
        return -get_log_likelihood (STRUC_DF, INDUC_DF, GAMMA, ALPHA, tau)
    
    start_time = time.time()
    if OPTION == 'basinhopping':
        alpha_max, gamma_max = optimize.basinhopping (ll, [0.5,0.5]).x
    elif OPTION == 'differential evolution':
        alpha_max, gamma_max = optimize.differential_evolution (ll, [(.1,0.99), (0.1,0.99)]).x
    elif OPTION == 'brute':
        alpha_max, gamma_max = optimize.brute (ll, [(0,1), (0,1)])[0]
    else:
        raise ValueError('Unknown option: {OPTION}')
    #print("--- %s seconds ---" % (time.time() - start_time))
    return alpha_max, gamma_max
