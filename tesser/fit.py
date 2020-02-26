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


def get_induction_log_likelihood(struc_df, induc_df, gamma, alpha, tau,
                                 return_trial=False, use_run=(2, 6)):
    """ This function gives the probability of obtaining the choices in the run, given
        specific values for alpha, gamma.
        INPUT:

        struc_df: Structured learning data in DataFrame format.
        induc_df: Generalized induction data in DataFrame format.
        gamma & alpha: discount and learning rate parameters. From 0.0 to 1.0.
    """

    # generate SR based on these parameters
    SR_all = sr.learn_sr(struc_df, gamma, alpha)
    SR = SR_all[use_run]

    # get likelihood of induction data
    num_trials = induc_df.shape[0]
    log_likelihood = 0
    all_trial_prob = np.zeros(num_trials)
    for i, trial in induc_df.iterrows():
        if np.isnan(trial.response):
            trial_probability = 0.5
            log_likelihood += np.log(trial_probability)
            all_trial_prob[i] = trial_probability
            continue
        trial_probability = probability_induction_choice(
            trial.cue, trial.opt1, trial.opt2, trial.response, SR, tau
        )
        eps = 0.000001
        if np.isnan(trial_probability):
            # probability undefined; can occur if SR has zeros
            trial_probability = eps
        elif trial_probability < eps:
            # probability too close to zero; set to minimal value
            trial_probability = eps
        log_likelihood += np.log(trial_probability)
        all_trial_prob[i] = trial_probability

    if return_trial:
        return log_likelihood, all_trial_prob
    else:
        return log_likelihood


def induction_brute (struc_df, induc_df):
    alphas = [float(i/20) for i in range(1, 20)]
    gammas = [float(i/20) for i in range(1, 20)]
    taus = [float(i/20) for i in range(1, 20)]
    likelihoods = [[[get_induction_log_likelihood(struc_df, induc_df, alphas[i], gammas[j], taus[k]) for k in range(19)] for j in range(19)]for i in range(19)]
    
    alpha_index, gamma_index, tau_index = 0, 0, 0
    likelihood_max = likelihoods[0][0]
    for i in range(19):
        for j in range(19):
            for k in range(19):
                if likelihoods[i][j][k] > likelihood_max:
                    alpha_index, gamma_index, tau_index = i, j, k
                    likelihood_max = likelihoods[i][j][k]
    
    return alphas[alpha_index], gammas[gamma_index], taus[tau_index]


def param_bounds(var_bounds, var_names):
    """Pack group-level parameters."""

    group_lb = [var_bounds[k][0] for k in var_names]
    group_ub = [var_bounds[k][1] for k in var_names]
    bounds = optimize.Bounds(group_lb, group_ub)
    return bounds


def function_logl_induct(struct_df, induct_df, fixed, var_names):
    """Generate log likelihood function for use with fitting."""

    param = fixed.copy()
    def fit_logl(x):
        param.update(dict(zip(var_names, x)))
        logl = get_induction_log_likelihood(struct_df, induct_df, **param)
        return -logl
    return fit_logl


def fit_induct(struct_df, induct_df, fixed, var_names, var_bounds,
               f_optim=optimize.differential_evolution,
               verbose=False, options=None):
    """Fit induction data."""

    if options is None:
        options = {}

    f_fit = function_logl_induct(struct_df, induct_df, fixed, var_names)
    bounds = param_bounds(var_bounds, var_names)
    res = f_optim(f_fit, bounds, disp=verbose, **options)

    # fitted parameters
    param = fixed.copy()
    param.update(dict(zip(var_names, res['x'])))

    logl = -res['fun']
    return param, logl


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
        
    elif option == 'induction brute':
        alpha_max, gamma_max, tau_max = induction_brute (struc_df, induc_df)
    else:
        raise ValueError('Unknown option: {option}')
    # print("--- %s seconds ---" % (time.time() - start_time))
    return alpha_max, gamma_max, tau_max


def grouping_error(struc_df, group_df, alpha, gamma):
    SR = sr.learn_sr(struc_df, gamma, alpha)[(1,5)]

    euclid_matrix = tasks.group_dist_mat(group_df)
    euclid_vector = distance.squareform(euclid_matrix)

    sr_matrix = util.make_sym_matrix(SR)
    sr_vector = distance.squareform(sr_matrix, checks=False) 

    slope, intercept, r_value, p_value, std_err = linregress(sr_vector, euclid_vector)
    
    euclid_estimates = np.array([slope*x + intercept for x in sr_vector])
    err = np.square(euclid_vector - euclid_estimates)
    err = np.mean(err)
    err = np.sqrt(err)
    
    return err

def group_brute (struc_df, group_df):
    alphas = [float(i/20) for i in range(1, 20)]
    gammas = [float(i/20) for i in range(1, 20)]
    errors = [[grouping_error(struc_df, group_df, alphas[i], gammas[j]) for j in range(19)] for i in range(19)]
    
    alpha_index, gamma_index = 0, 0
    error_min = errors[0][0]
    for i in range(19):
        for j in range(19):
            if errors[i][j] < error_min:
                alpha_index, gamma_index = i, j
                error_min = errors[i][j]
    
    return alphas[alpha_index], gammas[gamma_index]


def minimize_grouping_error(struc_df, group_df, option):
    def ge(x):
        alpha = x[0]
        gamma = x[1]
        return grouping_error(struc_df, group_df, alpha, gamma)

    if option == 'basinhopping':
        alpha_max, gamma_max = optimize.basinhopping(ge, [.5, .5]).x
    elif option == 'differential evolution':
        alpha_max, gamma_max = optimize.differential_evolution(ge, [(.01, 0.99), (0.01, 0.99)]).x
    elif option == 'brute':
        alpha_max, gamma_max = optimize.brute(ge, [(0.01, 1), (0.01, 1)])
    elif option == 'group brute':
        alpha_max, gamma_max = group_brute (struc_df, group_df)
    else:
        raise ValueError('Unknown option: {option}')
    # print("--- %s seconds ---" % (time.time() - start_time))
    return alpha_max, gamma_max
