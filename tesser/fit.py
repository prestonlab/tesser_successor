"""
Provides a set of functions designed to fit the parameters of the Successor
Representation learning agent to the induction and grouping data via maximum
likelihood.

pBGivenA computes the probability of choosing image B given a prompt image A,
using a learned SR matrix and Luce's choice rule.

pCGivenA computes the probability of the alternative choice.

probability_induction_choice combines these two functions to give a probability
of one of the agent's choices in the induction data.

get_induction_log_likelihood computes the joint probability of all of the subject's
choices present in the induction data given a learned SR matrix.

maximize_induction_likelihood finds parameters of the SR agent which maximize
the joint probability of all of the subject's choices under the SR matrix learned
by the agent.

"""


# model fitting/parameter optimization
import numpy as np
from . import util
from . import sr
from . import tasks
from scipy.spatial import distance
from scipy import optimize
from scipy.stats import linregress
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
    return SR[A][B] ** tau / (SR[A][B] ** tau + SR[A][C] ** tau)


def pCGivenA(A, B, C, SR, tau):
    """ Computes the probability of C given A using the Luce choice rule."""
    return 1 - pBGivenA(A, B, C, SR, tau)


def probability_induction_choice(cue, opt1, opt2, response, SR, tau):
    """ Computes the likelihood of a subject's response given the SR matrix for this subject."""
    if response == 0:
        return pBGivenA(cue, opt1, opt2, SR, tau)
    if response == 1:
        return pCGivenA(cue, opt1, opt2, SR, tau)


def get_induction_log_likelihood(struc_df, induc_df, gamma, alpha, tau, return_trial=False):
    """ This function gives the probability of obtaining the choices in the run, given
        specific values for alpha, gamma.
        INPUT:

        struc_df: Structured learning data in DataFrame format.
        induc_df: Generalized induction data in DataFrame format.
        gamma & alpha: discount and learning rate parameters. From 0.0 to 1.0.
    """
    SR = sr.explore_runs(struc_df, "once", gamma, alpha)
    SR_norm = SR / np.sum(SR)
    cue_sequence, opt1_sequence, opt2_sequence, response_sequence = util.get_induction_data(
        induc_df
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
        trial_probability = probability_induction_choice(
            cue_num, opt1_num, opt2_num, response_num, SR_norm, tau
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

    if return_trial:
        return log_likelihood, all_trial_prob
    else:
        return log_likelihood


def maximize_induction_likelihood(struc_df, induc_df, option):
    """ Numerically maximizes the log likelihood function on the set of the subject's choices
         to find optimal values for alpha, gamma
        INPUT:

        struc_df: Structured learning data in DataFrame format.
        induc_df: Generalized induction data in DataFrame format.
        option: scipy optimization functions:'basinhopping','differential evolution', 'brute'
    """

    def ll(x):
        alpha = x[0]
        gamma = x[1]
        tau = x[2]
        return -get_induction_log_likelihood(struc_df, induc_df, gamma, alpha, tau)

    start_time = time.time()
    if option == 'basinhopping':
        alpha_max, gamma_max, tau_max = optimize.basinhopping(ll, [0.5, 0.5, 0.5]).x
    elif option == 'differential evolution':
        alpha_max, gamma_max, tau_max = optimize.differential_evolution(ll,
                                                                        [(.000001, 0.99), (0.1, 0.99), (.00001, .99)]).x
    elif option == 'brute':
        alpha_max, gamma_max, tau_max = optimize.brute(ll, [(0, 1), (0, 1), (0, 1)])[0]
    else:
        raise ValueError('Unknown option: {option}')
    # print("--- %s seconds ---" % (time.time() - start_time))
    return alpha_max, gamma_max, tau_max


def grouping_error(struc_df, group_df, gamma, alpha):
    SR = sr.learn_sr(struc_df, gamma, alpha)
    SR = SR[2,6]
    euclid_matrix = np.array(group_df)
    euclid_vector = tasks.make_sym_matrix(euclid_matrix)
    sr_vector = tasks.make_sym_matrix(SR)
    slope, intercept, r_value, p_value, std_err = linregress(sr_vector[:, 0], euclid_vector[:, 0])
    return std_err


def minimize_grouping_error(struc_df, group_df, option):
    def ge(x):
        alpha = x[0]
        gamma = x[1]
        return grouping_error(struc_df, group_df, gamma, alpha)

    if option == 'basinhopping':
        alpha_max, gamma_max = optimize.basinhopping(ge, [.5, .5]).x
    elif option == 'differential evolution':
        alpha_max, gamma_max = optimize.differential_evolution(ge, [(.000001, 0.99), (0.1, 0.99)]).x
    elif option == 'brute':
        alpha_max, gamma_max = optimize.brute(ge, [(0, 1), (0, 1)])[0]
    else:
        raise ValueError('Unknown option: {option}')
    # print("--- %s seconds ---" % (time.time() - start_time))
    return alpha_max, gamma_max
