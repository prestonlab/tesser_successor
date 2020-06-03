import numpy as np
import emcee
import pymc3 as pm
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
    

def bayes_induct(struct_df, induct_df, num_samples, option):
    if option == 'emcee':
        fixed = {}
        var_names = ['gamma', 'alpha', 'tau']
        use_run = (1, 5)
        def log_prior(x):
            if min(x) < 0 or max(x) > 1:
                return -np.infty
            else:
                return 0
        def log_likelihood(x):
            return fit.get_induct_ll_all(struct_df, induct_df, fixed, var_names, x, use_run)
        def log_posterior(x):
            return log_prior(x) + log_likelihood(x)
        nwalkers = 50
        ndim = 3
        sampler = emcee.EnsembleSampler(nwalkers=nwalkers, ndim=ndim, log_prob_fn=log_posterior, args=[])
        x0 = uniform.rvs(size=(nwalkers, ndim))
        sampler.run_mcmc(x0, 10000)
        samples = sampler.get_chain(flat=True)
        return samples
    
    if option == 'pymc':
        fixed = {}
        var_names = ['gamma', 'alpha', 'tau']
        use_run = (1, 5)
        def log_likelihood(x):
            return fit.get_induct_ll_all(struct_df, induct_df, fixed, var_names, x, use_run)
        with pm.Model() as model:
            x = pm.Uniform('x', lower=0, upper=1, shape=3)
            y = pm.DensityDist('y', logp=log_likelihood, observed={'x': x})
            samples = pm.sample(10000)
        return samples
        
            
            
    