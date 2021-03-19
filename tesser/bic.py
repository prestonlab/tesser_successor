import numpy as np
import pandas as pd


def mean_stu_conf(array):
    mean = np.average(array)
    N =len(array)
    standard_error = np.std(array)/(N**0.5)
    return mean, standard_error


def params(results, error = False):
    d ={}
    for s in results.keys():
        if s == 'subject' or s == 'gamma2' or s == 'w' or s == 'w_prim' or s == 'w_bound1' or s == 'w_bound2':
            continue
        m, e =mean_stu_conf(results[s])
        d[s] =[m]
        if error:
            d[s[0]+'_err'] = [e]

    return d


def bic(k, n, L):
    ''' Function that creates a Bayesian information criterion .A model with the lowest BIC is preferred.
    n = number of data points
    k = number of free parameters
    L = Log liklihood
    '''
    return k*np.log(n) - 2*L
    
    
def get_bic(results, k, n, subject_BIC=False, array=True):
    # Sum of all log likelihoods 
    L = np.sum(results.log_like)
    
    subjects = results.subject
    # number of subjects, this creates the N multiplier for a group level BIC.
    N = len(subjects)
    
    # BIC at an individual level
    if subject_BIC:
        logs = np.array([results.query(f"subject == {sub}").log_like for sub in subjects])
        BIC = [bic(k,n, l[0]) for l in logs]
        if array:
            bics = np.array(BIC)
        else:
            bics = pd.DataFrame({'BIC':BIC},subjects)
        return bics
    else:
        # BIC at a group level
        return bic(N*k, N*n, L)
    

def wbic(a, axis=1):
    """Bayesian weights."""
    min_bic = np.expand_dims(np.min(a, axis), axis)
    delta_bic = np.exp(-0.5 * (a - min_bic))
    sum_bic = np.expand_dims(np.sum(delta_bic, axis), axis)
    return delta_bic / sum_bic
