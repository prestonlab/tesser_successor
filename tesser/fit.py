import numpy as np
import sr

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
