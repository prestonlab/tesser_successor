"""
Provides a set of functions designed to fit the parameters of the Successor
Representation learning agent to the induction and grouping data via maximum
likelihood.

pBGivenA computes the probability of choosing image B given a prompt image A,
using a learned SR matrix and Luce's choice rule.

pCGivenA computes the probability of the alternative choice.

probability_induction_choice combines these two functions to give a probability
of one of the agent's choices in the induction data.

get_induction_log_likelihood computes the joint probability of all of the
subject's choices present in the induction data given a learned SR matrix.

maximize_induction_likelihood finds parameters of the SR agent which maximize
the joint probability of all of the subject's choices under the SR matrix
learned by the agent.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from . import util
from . import sr
from . import tasks
from . import network
from scipy.spatial import distance
from scipy import optimize
from scipy.stats import linregress
import time


def eu_dist(a, b, SR):
    return distance.euclidean(SR[a], SR[b])

def prob_induct_choice_hybrid(cue, opt, response, SR, comm, w, tau, choice_rule):
    """Likelihood of induction response."""

    if np.all(SR[cue, opt] == 0):
        return np.nan

    support1 = w * SR[cue, opt[0]] + (1 - w) * comm[cue, opt[0]]
    support2 = w * SR[cue, opt[1]] + (1 - w) * comm[cue, opt[1]]
    support = [support1, support2]
    if choice_rule == 'softmax':
        prob = (np.exp(support[response] / tau) /
               (np.exp(support[0] / tau) + np.exp(support[1] / tau)))
    else:
        prob = ((support[response] ** tau) /
            (support[0] ** tau + support[1] ** tau))
    return prob


def prob_induct_choice(cue, opt, response, SR, tau):
    """Likelihood of induction response."""

    if np.all(SR[cue, opt] == 0):
        return np.nan

    prob = ((SR[cue, opt[response]] ** tau) /
            (SR[cue, opt[0]] ** tau + SR[cue, opt[1]] ** tau))
    return prob


def prob_induct_subject(struct, induct, gamma, alpha, tau, w=1,
                        response_key='response', use_run=(2, 6),comm= [], choice_rule='softmax'):
    """Calculate induction task probabilities for one subject."""

    # generate SR based on these parameters
    SR_all = sr.learn_sr(struct, gamma, alpha)
    SR = SR_all[use_run]
    induct = induct.reset_index()



    # get likelihood of induction data
    num_trials = induct.shape[0]
    trial_prob = np.zeros(num_trials)
    for i, trial in induct.iterrows():
        if np.isnan(trial.response):
            trial_prob[i] = np.nan
            continue

        trial_prob[i] = prob_induct_choice_hybrid(
            trial.cue, [trial.opt1, trial.opt2], int(trial[response_key]), SR, comm, w, tau, choice_rule)
    return trial_prob


def assess_induct_fit_subject(struct, induct, param):
    """Compare model and data in fitting the induction task."""

    trial_prob = prob_induct_subject(struct, induct, param['gamma'],
                                     param['alpha'], param['tau'],
                                     response_key='response')
    induct = induct.copy()
    induct.loc[:, 'Data'] = induct['Acc']
    induct.loc[:, 'Model'] = trial_prob
    results = induct.melt(id_vars=['SubjNum', 'TrialNum', 'QuestType',
                                   'Environment'], value_vars=['Data', 'Model'],
                          var_name='Source', value_name='Accuracy')
    return results


def assess_induct_fit_subject_hybrid(struct, induct, param, comm, choice_rule):
    """Compare model and data in fitting the induction task."""

    trial_prob = prob_induct_subject(struct, induct, param['gamma'],
                                     param['alpha'], param['tau'], param['w'],
                                     response_key='response', comm=comm, choice_rule=choice_rule)
    induct = induct.copy()
    induct.loc[:, 'Data'] = induct['Acc']
    induct.loc[:, 'Model'] = trial_prob
    results = induct.melt(id_vars=['SubjNum', 'TrialNum', 'QuestType',
                                   'Environment'], value_vars=['Data', 'Model'],
                          var_name='Source', value_name='Accuracy')
    return results


def plot_induct_fit(results):
    """Plot fit to induction performance."""

    fig, ax = plt.subplots(1, 2)
    sns.pointplot(x='QuestType', y='Accuracy', hue='Source',
                  data=results, dodge=True, ax=ax[0])
    sns.pointplot(x='Environment', y='Accuracy', hue='Source',
                  data=results, dodge=True, ax=ax[1])


def get_induction_log_likelihood_hybrid(struc_df, induc_df, gamma, alpha, tau, w = 1.0,
                                 return_trial=False, use_run=(2, 6), comm=[], choice_rule='softmax'):
    """ This function gives the probability of obtaining the choices in the run,
        given specific values for alpha, gamma.
        INPUT:

        struc_df: Structured learning data in DataFrame format.
        induc_df: Generalized induction data in DataFrame format.
        gamma & alpha: discount and learning rate parameters. From 0.0 to 1.0.
        tau is a waiting parameter in SR. Larger than 0.0.
        w is a waiting parameter for SR in a different community matrix. From 0.0 to 1.0.
            0 is for a different community matrix
            1 is for the same community matrix
    """

    # generate SR based on these parameters
    SR_all = sr.learn_sr(struc_df, gamma, alpha)
    SR = SR_all[use_run]

    induc_df = induc_df.reset_index()
    if comm == []:
        # get community matrix
        net = network.temp_node_info()
        comm = 1 - distance.squareform(distance.pdist(net['comm'][:, None], 'hamming'))
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
        trial_probability = prob_induct_choice_hybrid(
            trial.cue, [trial.opt1, trial.opt2], int(trial.response), SR, comm, w, tau, choice_rule)
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

    
def get_induct_ll_all(struct_df, induct_df, fixed, var_names, x, use_run, comm=[], choice_rule='softmax'):
    param = fixed.copy()
    flexible = {}
    flex_names = []
    for var_name in var_names:
        if var_name not in list(fixed.keys()):
            flex_names.append(var_name)

    num_flex = len(flex_names)
    for i in range(num_flex):
        key = flex_names[i]
        flexible[key] = x[i]
    param.update(flexible)
    logl = 0
    subjects = struct_df.SubjNum.unique()
    for subject in subjects:
        subj_filter = f'SubjNum == {subject}'
        subj_struct = struct_df.query(subj_filter)
        subj_induct = induct_df.query(subj_filter)
        subj_logl = get_induction_log_likelihood_hybrid(subj_struct,subj_induct, **param,
                                                 return_trial=False, use_run=use_run, comm=comm, choice_rule=choice_rule)
        logl += subj_logl
    
    return logl


def induction_brute(struc_df, induc_df):
    alphas = [float(i / 20) for i in range(1, 20)]
    gammas = [float(i / 20) for i in range(1, 20)]
    taus = [float(i / 20) for i in range(1, 20)]
    likelihoods = [[[get_induction_log_likelihood_hybrid(
        struc_df, induc_df, alphas[i], gammas[j], taus[k])
        for k in range(19)] for j in range(19)] for i in range(19)]

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


def fit_induct(struct_df, induct_df, fixed, var_names, var_bounds,
               f_optim=optimize.differential_evolution,
               verbose=False, options=None, use_run=(2, 6), comm=[], choice_rule='softmax'):
    """Fit induction data for one subject.

    For a given set of parameters, the structure learning task is used
    to generate a simulated SR matrix. Then this matrix is used to 
    simulate responses in the induction task. Parameters are optimized
    to obtain the set that maximizes the probability of the responses
    observed in the induction task.

    Parameters
    ----------
    struct_df : DataFrame
        Structure learning data.

    induct_df : DataFrame
        Induction test data.

    fixed : dict
        Parameter values for all fixed parameters.

    var_names : list
        String name for each variable parameter.

    var_bounds : dict
        Bounds (in low, high order) for each variable parameter.

    f_optim : function
        Function to use for parameter optimization.

    verbose : Boolean
        If true, more information about the search will be printed.

    options : dict
        Options to pass to f_optim.

    use_run : tuple
        Run to take the SR matrix from for predicting induction data,
        specified as (part_number, run_number).

    Returns
    -------
    param : dict
        Best-fitting parameters.

    logl : float
        Maximum log likelihood.
    """

    if options is None:
        options = {}

    param = fixed.copy()
    subjects = struct_df.SubjNum.unique()

    def f_fit(x):
        param.update(dict(zip(var_names, x)))
        # draft code to fit at the group level:
        logl = 0
        for subject in subjects:
            subj_filter = f'SubjNum == {subject}'
            subj_struct = struct_df.query(subj_filter)
            subj_induct = induct_df.query(subj_filter)
            subj_logl = get_induction_log_likelihood_hybrid(subj_struct,subj_induct, **param,
                                                     return_trial=False, use_run=use_run, comm=comm, choice_rule=choice_rule)
            logl += subj_logl
        # logl = get_induction_log_likelihood(struct_df, induct_df, **param,
        #                                     return_trial=False, use_run=use_run)
        return -logl

    bounds = param_bounds(var_bounds, var_names)
    res = f_optim(f_fit, bounds, disp=verbose, **options)

    # fitted parameters
    param = fixed.copy()
    param.update(dict(zip(var_names, res['x'])))

    logl = -res['fun']
    return param, logl


def fit_induct_indiv(struct, induct, fixed, var_names, var_bounds, comm=[], choice_rule='softmax'):
    """Estimate parameters for individual subjects."""

    df_list = []
    for sub in struct['SubjNum'].unique():
        print(f'Estimating parameters for {sub}...')
        subj_struct = struct.query(f'SubjNum == {sub}')
        subj_induct = induct.query(f'SubjNum == {sub}')
        param, logl = fit_induct(subj_struct, subj_induct, fixed,
                                 var_names, var_bounds, verbose=False,comm=comm, choice_rule=choice_rule)
        param['subject'] = sub
        param['log_like'] = logl
        df = pd.DataFrame(param, index=[0])
        df_list.append(df)
    df = pd.concat(df_list, axis=0, ignore_index=True)
    #df = df.set_index('subject')
    return df


def maximize_induction_likelihood(struc_df, induc_df, option):
    """ Numerically maximizes the log likelihood function on the set of
        the subject's choices to find optimal values for alpha, gamma
        INPUT:

        struc_df: Structured learning data in DataFrame format.
        induc_df: Generalized induction data in DataFrame format.
        option: scipy optimization functions:'basinhopping',
        'differential evolution', 'brute'
    """

    def ll(x):
        alpha = x[0]
        gamma = x[1]
        tau = x[2]
        return -get_induction_log_likelihood(struc_df, induc_df, gamma,
                                             alpha, tau)

    start_time = time.time()
    if option == 'basinhopping':
        alpha_max, gamma_max, tau_max = optimize.basinhopping(
            ll, [0.5, 0.5, 0.5]).x
    elif option == 'differential evolution':
        alpha_max, gamma_max, tau_max = optimize.differential_evolution(
            ll, [(.000001, 0.99), (0.1, 0.99), (.00001, .99)]).x
    elif option == 'brute':
        alpha_max, gamma_max, tau_max = optimize.brute(
            ll, [(0, 1), (0, 1), (0, 1)])[0]

    elif option == 'induction brute':
        alpha_max, gamma_max, tau_max = induction_brute(struc_df, induc_df)
    else:
        raise ValueError('Unknown option: {option}')
    # print("--- %s seconds ---" % (time.time() - start_time))
    return alpha_max, gamma_max, tau_max


def grouping_error(struc_df, group_df, alpha, gamma):
    SR = sr.learn_sr(struc_df, gamma, alpha)[(1, 5)]

    euclid_matrix = tasks.group_dist_mat(group_df)
    euclid_vector = distance.squareform(euclid_matrix)

    sr_matrix = util.make_sym_matrix(SR)
    sr_vector = distance.squareform(sr_matrix, checks=False)

    slope, intercept, r_value, p_value, std_err = linregress(sr_vector,
                                                             euclid_vector)

    euclid_estimates = np.array([slope * x + intercept for x in sr_vector])
    err = np.square(euclid_vector - euclid_estimates)
    err = np.mean(err)
    err = np.sqrt(err)

    return err


def group_brute(struc_df, group_df):
    alphas = [float(i / 20) for i in range(1, 20)]
    gammas = [float(i / 20) for i in range(1, 20)]
    errors = [[grouping_error(struc_df, group_df, alphas[i], gammas[j])
               for j in range(19)] for i in range(19)]

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
        alpha_max, gamma_max = optimize.differential_evolution(
            ge, [(.01, 0.99), (0.01, 0.99)]).x
    elif option == 'brute':
        alpha_max, gamma_max = optimize.brute(ge, [(0.01, 1), (0.01, 1)])
    elif option == 'group brute':
        alpha_max, gamma_max = group_brute(struc_df, group_df)
    else:
        raise ValueError('Unknown option: {option}')
    # print("--- %s seconds ---" % (time.time() - start_time))
    return alpha_max, gamma_max

