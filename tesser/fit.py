import numpy as np
import sr

def pBGivenA (A, B, C, SR):
    return SR[A][B] / (SR[A][B] + SR[A][C])

def pCGivenA (A, B, C, SR):
    return 1 - pBGivenA (A, B, C, SR)

def likelihood (cue, opt1, opt2, response, SR):
    if response == 0:
        return pBGivenA (cue, opt1, opt2, SR)
    if response == 1:
        return pCGivenA (cue, opt1, opt2, SR)
    

def get_log_likelihood(SUBJECT, GAMMA, ALPHA):
    ''' This function gives the probability of obtaining the choices in the run, given
    specific values for alpha, gamma.
    '''
    SR = sr.explore_runs ('once', SUBJECT, GAMMA, ALPHA)
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

def maximize_likelihood (subject):
    h = 10e-3
    alpha_max, gamma_max = 0.0, 0.0
    alpha, gamma = h, h
    log_likelihood_max = get_log_likelihood (subject, gamma, alpha)
    log_likelihoods = []
    
    
    while alpha <= 1:
        gamma = h
        while gamma <= 1:
            log_likelihood = get_log_likelihood(alpha, gamma, obj_sequence_directory, induction_directory)
            log_likelihoods.append (log_likelihood)
            if log_likelihood > log_likelihood_max:
                log_likelihood_max = log_likelihood
                alpha_max, gamma_max = alpha, gamma
            gamma += h
        alpha += h
    
    return alpha_max, gamma_max
