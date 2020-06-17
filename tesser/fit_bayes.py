import numpy as np
import emcee
import pymc3 as pm
import pystan as ps
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
    
    if option == 'stan':
        struct_arr = struct_df.objnum
        cue = induct_df.CueNum.to_numpy()
        response1 = induct_df.Opt1Num.to_numpy()
        response2 = induct_df.Opt2Num.to_numpy()
        choice = induct_df.response.to_numpy().astype(int)
        indices = np.logical_or(choice == 1, choice == 0)
        cue = cue[indices]
        response1 = response1[indices]
        response2 = response2[indices]
        choice = choice[indices]
        num_struct = struct_arr.shape[0]
        num_induct = cue.shape[0]
        model = """
        data {
            int<lower=0> num_struct;
            int<lower=0> num_induct;
            int struct_arr[num_struct];
            int cue[num_induct];
            int response1[num_induct];
            int response2[num_induct];
            int choice[num_induct];
        }
        parameters {
            real<lower=0, upper=1> alpha;
            real<lower=0, upper=1> gam;
            real<lower=0> tau;
        }
        model {
            matrix[21, 21] SR = rep_matrix(0, 21, 21);
            vector[21] ones = rep_vector(1, 21);
            matrix[21, 21] I = diag_matrix(ones);
            int s;
            int s_new;
            int A;
            int B;
            int C;
            real p;
            alpha ~ uniform(0, 1);
            gam ~ uniform(0, 1);
            tau ~ gamma(2, 4);
            for (i in 1:num_struct-1){
                s = struct_arr[i];
                s_new = struct_arr[i+1];
                SR[s] = (1 - alpha) * SR[s] + alpha * (
                I[s_new] + gam * SR[s_new]);
            }
            for (i in 1:num_induct){
                A = cue[i];
                B = response1[i];
                C = response2[i];
                p = exp(SR[A][B]/tau) / (exp(SR[A][B]/tau) + exp(SR[A][C]/tau));
                choice[i] ~ bernoulli(p);
            }
        }
        """
        data = {'num_struct': num_struct, 'num_induct': num_induct, 'struct_arr': struct_arr, 'cue': cue, 'response1': response1, 
               'response2': response2, 'choice': choice}
        sm = ps.StanModel(model_code=model)
        samples = sm.sampling(data=data, iter=num_samples, chains=4, warmup=500, thin=1, seed=101)
        return samples
    
    
            
    