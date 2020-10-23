cimport cython
from libc.math cimport exp
from libc.math cimport log
from libc.math cimport isnan
import numpy as np


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True) 
cpdef prob_choice(int cue, 
                   int opt1, 
                   int opt2, 
                   int response, 
                   double [:,:] SR,
                   double [:,:] model,
                   double w,
                   double tau):
    
    cdef double support1
    cdef double support2
    cdef double support
    cdef double prob


    
    support1 = w * SR[cue, opt1] + (1.0 - w) * model[cue, opt1]
    support2 = w * SR[cue, opt2] + (1.0 - w) * model[cue, opt2]
    if response == 0:
        support = support1
    else:
        support = support2
    prob = (exp(support/ tau) / (exp(support1 / tau) + exp(support2 / tau)))

    return prob


@cython.boundscheck(False)
@cython.wraparound(False)
def prob_induct(int [:] cue,
                  int [:] opt1,
                  int [:] opt2,
                  int[:] response,
                  double [:,:] SR,
                  double [:,:] model,
                  double w,
                  double tau,
                  bint return_trial,
                  double [:] all_trial_prob):
    
    
    cdef Py_ssize_t num_trials = cue.shape[0]
    cdef double log_likelihood = 0.0
    cdef double trial_probability
    cdef double eps
    
        
    for i in range(num_trials):
        if np.isnan(response[i]):
            trial_probability = 0.5
            log_likelihood += log(trial_probability)
            all_trial_prob[i] = trial_probability
            continue
        trial_probability = prob_choice(
            cue[i], opt1[i], opt2[i], response[i], SR, model, w, tau)
        eps = 0.000001
        if np.isnan(trial_probability):
            # probability undefined; can occur if SR has zeros
            trial_probability = eps
        elif trial_probability < eps:
            # probability too close to zero; set to minimal value
            trial_probability = eps
        log_likelihood += log(trial_probability)
        all_trial_prob[i] = trial_probability

    if return_trial:
        return log_likelihood, all_trial_prob
    else:
        return log_likelihood


@cython.boundscheck(False)
@cython.wraparound(False)
def prob_induct_subject(double [:,:] SR,
                  int [:] cue,
                  int [:] opt1,
                  int [:] opt2,
                  int[:] response,
                  double tau,
                  double w,
                  double [:,:] model,
                  double [:] trial_prob):
    
    
    cdef Py_ssize_t num_trials = cue.shape[0]    
        
    for i in range(num_trials):
        if np.isnan(response[i]):
            trial_prob[i] = np.nan
            continue
            
        trial_prob[i] = prob_choice(
            cue[i], opt1[i], opt2[i], response[i], SR, model, w, tau)
    return trial_prob
