import numpy as np
from tesser import fit
from scipy.stats import multivariate_normal as mvn
from scipy.stats import uniform


class MHSampler:
    def __init__(self, f, num_burn, num_save, d, x0):
        self.f = f
        self.num_burn = num_burn
        self.num_save = num_save
        self.d = d
        self.samples = np.zeros((num_burn+num_save, d))
        self.samples[0] = x0
    
    def sample(self):
        n = self.num_burn + self.num_save
        for t in range(n-1):
            x_new = mvn.rvs(mean=self.samples[t], cov=.25)
            log_acceptance = self.f(x_new) - self.f(self.samples[t])
            log_correction = np.log(mvn.pdf(self.samples[t], mean=x_new, cov=.25))
            log_correction -= np.log(mvn.pdf(x_new, mean=self.samples[t], cov=.25))
            log_acceptance += log_correction
            log_acceptance = min(0, log_acceptance)
            u = uniform.rvs()
            if np.log(u) <= log_acceptance:
                self.samples[t+1] = x_new
            else:
                self.samples[t+1] = self.samples[t]
            if t % 10 == 0:
                print ("Status: t >", t)
    
    def get_samples(self):
        return self.samples[self.num_burn:]
    

def bayes_induct(struct_df, induct_df, fixed, var_names, prior, use_run, num_burn, num_save):
    d = len(var_names) - len(list(fixed.keys()))
    def f(x):
        if prior(x) == 0:
            return -np.infty
        else:
            return fit.get_induct_ll_all(struct_df, induct_df, fixed, var_names, x, use_run) + np.log(prior(x))
    mhSampler = MHSampler(f, num_burn, num_save, d, np.array([.5 for i in range(d)]))
    mhSampler.sample()
    samples = mhSampler.get_samples()
    return samples
    
    