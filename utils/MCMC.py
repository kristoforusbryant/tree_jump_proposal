import numpy as np
import scipy.stats as st
from tqdm import tqdm
import copy

def MCMC(prior_s, prior_r_e, lik_r_e, prop_s, prop_r_e, data, N=5000, burnin=2500): 
    """
    prior_s(f): Samples from prior distribution 
    prior_e(f): Evaluates prior probability of current sample 
    lik_e(f): Evaluates likelihood probability of current sample
    prop_s(f): Samples from proposal distribution given the previous sample  
    prop_e(f): Evaluates probability of the current sample given the previous sample 
    """
    
    results = {'SAMPLES':[],
               'ALPHAS':[],
               'PARAMS':[], 
               'INDEX':[], 
               'LIK_R':[], 
               'PRIOR_R':[], 
               'PROP_R':[]
              }
    
    # Initialisation 
    params = prior_s()
    print(params)
    
    # Burn in 
    for i in tqdm(range(burnin)):
        results['PARAMS'].append(params)
        params_ = prop_s(params)
        
        lik_ratio = lik_r_e(data, params, params_)
        prior_ratio = prior_r_e(params, params_)
        prop_ratio = prop_r_e(params, params_)
        
        results['LIK_R'].append(lik_ratio) 
        results['PRIOR_R'].append(prior_ratio)
        results['PROP_R'].append(prop_ratio) 
        
        alpha = min(lik_ratio + prior_ratio + prop_ratio, 0)
        
        results['ALPHAS'].append(alpha)
        if np.log(np.random.uniform()) < alpha: 
            params = params_ 
            params_copy = copy.deepcopy(params)
            
    # Sampler
    for i in tqdm(range(N)):
        results['SAMPLES'].append(params)
        params_ = prop_s(params)
        results['PARAMS'].append(params_)
        
        lik_ratio = lik_r_e(data, params, params_)
        prior_ratio = prior_r_e(params, params_)
        prop_ratio = prop_r_e(params, params_)
        
        results['LIK_R'].append(lik_ratio) 
        results['PRIOR_R'].append(prior_ratio)
        results['PROP_R'].append(prop_ratio) 
        
        alpha = min(lik_ratio + prior_ratio + prop_ratio, 0)
        
        if np.isnan(alpha):
            print("!!!! ALPHA IS NAN !!!!")
            print(lik_r_e(data, params_))
            print(params_)
            
        results['ALPHAS'].append(alpha)
        if np.log(np.random.uniform()) < alpha: 
            params = copy.deepcopy(params_) 
            
    return results

class Sampler(): 
    def __init__(self, param_l, fun_l):
        assert len(param_l) == len(fun_l), "number of parameter-sets and functions do not match"
        self.param_l = param_l 
        self.fun_l = fun_l 
        
    def __call__(self): 
        sample = []
        for i in range(len(self.param_l)): 
            sample.append(self.fun_l[i](*self.param_l[i]))
        return tuple(sample)
    
class Evaluator(): 
    def __init__(self, param_l, fun_l):
        assert len(param_l) == len(fun_l), "number of parameter-sets and functions do not match"
        self.param_l = param_l 
        self.fun_l = fun_l 
        
    def __call__(self,  data): 
        assert len(data) == len(self.fun_l), "number of data and functions do not match"
        maxmin = lambda x: max(min(x,np.exp(700)),np.exp(-700))
        sample = []
        for i in range(len(self.param_l)):
            sample.append(np.log(maxmin(self.fun_l[i](data[i], *self.param_l[i]))))
        return tuple(sample)
