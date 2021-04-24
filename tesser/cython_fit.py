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
from . import cython_sr
from . import sr
from . import tasks
from . import network
from . import cfit
from . import csr
from scipy.spatial import distance
from scipy import optimize
from scipy.stats import linregress



def plot_induct_fit(results):
    """Plot fit to induction performance."""

    fig, ax = plt.subplots(1, 2)
    sns.pointplot(x='QuestType', y='Accuracy', hue='Source',
                  data=results, dodge=True, ax=ax[0])
    sns.pointplot(x='Environment', y='Accuracy', hue='Source',
                  data=results, dodge=True, ax=ax[1])

    
def param_bounds(var_bounds, var_names):
    """Pack group-level parameters."""

    group_lb = [var_bounds[k][0] for k in var_names]
    group_ub = [var_bounds[k][1] for k in var_names]
    bounds = optimize.Bounds(group_lb, group_ub)
    return bounds



def assess_induct_fit_subject_hybrid(struct, induct, param, n_states, split, model_type, model=[]):
    """Compare model and data in fitting the induction task."""
    subject = int(struct.SubjNum.unique())
    SR = cython_sr.learn_sr(struct, param["gamma"], param["alpha"], n_states)
    if model_type =='comm':
        model_subj = model
    if (model_type == 'multiple gamma') or (model_type == 'multiple zero'):
        model_subj = cython_sr.learn_sr(struct, param["gamma2"], param["alpha"], n_states)
    if model_type =='true transitional':
        model_subj = model[subject]
    else:
        model_subj = np.zeros([n_states,n_states])
        

    induct_df = induct.reset_index()
    num_trials = induct_df.shape[0]
    cue = induct_df.cue.to_numpy()
    cue = cue.astype(np.dtype('i'))
    opt1 = induct_df.opt1.to_numpy()
    opt1 = opt1.astype(np.dtype('i'))
    opt2 = induct_df.opt2.to_numpy()
    opt2 = opt2.astype(np.dtype('i'))
    res = induct_df.response.to_numpy()
    res = res.astype(np.dtype('i'))
    trial_prob = np.zeros(num_trials)
    questions = induct.QuestType.unique()
    if split:
        trial_prob = 0
        for question in questions:
            for p in param.keys():
                if question.lower() in p:
                    
                    w = param[p]
                    trial_prob += cfit.prob_induct_subject(SR, cue, opt1, opt2, res ,
                                      param['tau'], w,
                                      model=model_subj, trial_prob=trial_prob)
            
    else:
#         print(type(cue)) # use to debug 
        trial_prob = cfit.prob_induct_subject(SR, cue, opt1, opt2, res ,
                                          param['tau'], param['w'],
                                          model=model_subj, trial_prob=trial_prob)
    induct = induct_df.copy()
    induct.loc[:, 'Data'] = induct['Acc']
    induct.loc[:, 'Model'] = trial_prob
    results = induct.melt(id_vars=['SubjNum', 'TrialNum', 'QuestType',
                                   'Environment'], value_vars=['Data', 'Model'],
                          var_name='Source', value_name='Accuracy')
    return results


def get_induction_log_likelihood_hybrid(SR, induc_df, tau, w, n_states, return_trial, model):
    

    induc_df = induc_df.reset_index()
    num_trials = induc_df.shape[0]
    cue = induc_df.cue.to_numpy()
    cue = cue.astype(np.dtype('i'))
    opt1 = induc_df.opt1.to_numpy()
    opt1 = opt1.astype(np.dtype('i'))
    opt2 = induc_df.opt2.to_numpy()
    opt2 = opt2.astype(np.dtype('i'))
    res = induc_df.response.to_numpy()
    res = res.astype(np.dtype('i'))
    all_trial_prob = np.zeros(num_trials)
    return cfit.prob_induct(cue, opt1, opt2, res , SR, model, w, tau, return_trial, all_trial_prob)

    
def fit_induct(struct_df, induct_df, fixed, split, var_names, var_bounds, n_states,
               f_optim=optimize.differential_evolution,
               verbose=False, options=None, model_type='comm', model=[], split_list = []):
    """Fit induction data for one subject faster cython improved version.

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
    questions = induct_df.QuestType.unique()

    def f_fit(x):
        param.update(dict(zip(var_names, x)))
        # draft code to fit at the group level:
        logl = 0
        for subject in subjects:
            subj_filter = f'SubjNum == {subject}'
            subj_struct = struct_df.query(subj_filter)
            subj_induct = induct_df.query(subj_filter)
            
            SR = cython_sr.learn_sr(subj_struct, param["gamma"], param["alpha"], n_states)
            if model_type == 'multiple gamma':
                model_subj = cython_sr.learn_sr(subj_struct, param["gamma2"], param["alpha"], n_states)

            if model_type =='true transitional':
                model_subj = model[subject]
                
            if model_type =='comm':
                model_subj = model
                
            else:
                model_subj = np.zeros([n_states,n_states])
            
            if split:
                subj_logl = 0
                for question in questions:
                    question_param = param.copy()
                    question_induct = subj_induct.query(f'QuestType == "{question}"')
                    for sp in split_list:
                        param_name = f"{sp}_{question.lower()}"

                        question_param[sp] = param[param_name]
                        
                    w = question_param["w"]
                    tau = question_param["tau"]
                    subj_logl += get_induction_log_likelihood_hybrid(SR, question_induct, tau, w, n_states,
                                                     False, model_subj)
            else:
                w = param["w"]
                tau = param["tau"]
                subj_logl = get_induction_log_likelihood_hybrid(SR, subj_induct, tau, w, n_states,
                                                     False, model_subj)

            logl += subj_logl
        # logl = get_induction_log_likelihood_hybrid(SR, subj_induct, tau = param["tau"], param["w"], n_states=n_states,
#                                                          return_trial=False, model=model_subj)
        return -logl

    bounds = param_bounds(var_bounds, var_names)
    res = f_optim(f_fit, bounds, disp=verbose, **options)

    # fitted parameters
    param = fixed.copy()
    param.update(dict(zip(var_names, res['x'])))

    logl = -res['fun']
    return param, logl


def fit_induct_indiv(struct, induct, fixed, var_names, var_bounds, split, n_states, verbose, model_type, model, split_list, free_params):
    """Estimate parameters for individual subjects."""

    df_list = []
    for sub in struct['SubjNum'].unique():
        print(f'Estimating parameters for {sub}...')
        subj_struct = struct.query(f'SubjNum == {sub}')
        subj_induct = induct.query(f'SubjNum == {sub}')
        param, logl = fit_induct(
                subj_struct, subj_induct, fixed, split, var_names, var_bounds, n_states,
                f_optim=optimize.differential_evolution, verbose=False, options=None,
                model_type=model_type, model=model, split_list=split_list)
        param['subject'] = sub
        param['log_like'] = logl
        df = pd.DataFrame(param, index=[0])
        df_list.append(df)
    df = pd.concat(df_list, axis=0, ignore_index=True)
    df['k'] = free_params
    #df = df.set_index('subject')
    return df

def fitted_results(struct_all, induct_all, results, n_states, split, model_type, model=[]):
    ''' Results for accuracy of individual fitting by question type.
        INPUT:
            results: DataFrame with param
        OUTPUT:
            fitted: Results fitted by accuracy
    '''
    res_list = []
    for subject in results.index.unique():
        subj_filter = f'SubjNum == {subject}'
        subj_struct = struct_all.query(subj_filter)
        subj_induct = induct_all.query(subj_filter)
        subj_param = results.loc[subject]

        param={}
        for name in results.columns:
            if (name == 'log_like') or (name == 'k'):
                continue

            param[name] = subj_param[name]
        res = assess_induct_fit_subject_hybrid(subj_struct , subj_induct, param, n_states, split, model_type, model)
        res_list.append(res)
    fitted = pd.concat(res_list, axis=0)
    return fitted


def plot_by_question(fitted, split, fig_name):
    '''Plots fitted results for accuracy of individual fitting by question type.
        INPUT:
            fitted: DataFrame with paramarmets fitted by accuracy
        OUTPUT:
            fig: Accuracy plots by question type
    '''
    ########################################################################### make into function
#     res_list = []
#     for subject in results.index.unique():
#         subj_filter = f'SubjNum == {subject}'
#         subj_struct = struct_all.query(subj_filter)
#         subj_induct = induct_all.query(subj_filter)
#         subj_param = results.loc[subject]

# #         print(subject)
#         param={}
#         for name in results.columns:
#             if name == 'log_like':
#                 continue

#             param[name] = subj_param[name]
#         res = assess_induct_fit_subject_hybrid(subj_struct , subj_induct, param, n_states, split, model_type, model)
#         res_list.append(res)
#     fitted = pd.concat(res_list, axis=0)
    
    ###########################################################################
    fig_name = fig_name.upper()
    if split:
        fig_name += '_weighted' 
    fig, axes = plt.subplots(2, 2, figsize = (10,10),sharey=False)
    fig.suptitle(fig_name) 
    names = [ 'Environment', 'QuestType']
    for i, name in enumerate(names):
        n = fitted.groupby(['Source', 'SubjNum', name])['Accuracy'].mean()
        m  = n.unstack(level=0)
        g = sns.scatterplot(x='Model', y='Data', hue=name, data=m.reset_index(), ax=axes[0][i % 2])
        g.set_xlim(0, 1.02)
        g.set_ylim(0, 1.02)
        g.set_aspect(1)
        g.set_title('Individual accuracies \n data vs model \n by '+name)
        g.plot((0, 1), (0, 1), '-k')

        f = sns.pointplot(kind='point', x=name, y='Accuracy', 
                    hue='Source', dodge=True, data=n.reset_index(), ax=axes[1][i % 2])
        f.set(ylim=(0, 1.02))
#     fig.savefig(path+fig_name)
    return fig

############################################################
########### python setup.py build_ext --inplace ############
############################################################
