import numpy as np
import scipy.stats as st
from tqdm import tqdm
import copy

def MCMC(prior_s, prior_e, lik_e, prop_s, prop_e, data, N=1000, burnin=500): 
    """
    prior_s(f): Samples from prior distribution 
    prior_e(f): Evaluates prior probability of current sample 
    lik_e(f): Evaluates likelihood probability of current sample
    prop_s(f): Samples from proposal distribution given the previous sample  
    prop_e(f): Evaluates probability of the current sample given the previous sample 
    """
    
    l = []
    alphas = []
    params_l = []
    post_l = []
    idx = []
    
    maxmin = lambda x: max(min(x,700),-700)
    post_e = lambda lik,prior,params: lik(data, params) + np.sum(prior(params))
    # Initialisation 
    while True: 
        params = prior_s()
        if post_e(lik_e, prior_e, params) != -np.inf: 
            break
    
    # Burn in 
    for i in tqdm(range(burnin)):
        params_l.append(params)
        params_ = prop_s(params)
        numerator = post_e(lik_e, prior_e, params_) + prop_e(params,params_) 
        denominator = post_e(lik_e, prior_e, params) + prop_e(params_,params)
        post_l.append((numerator, denominator)) 
        alpha = min(numerator - denominator, 0)
        
        alphas.append(alpha)
        if np.random.uniform() < np.exp(maxmin(alpha)): 
            params = params_ 
            params_copy = copy.deepcopy(params)
            
    # Sampler
    for i in tqdm(range(N)):
        params_l.append(params)
        params_ = prop_s(params)
        numerator = post_e(lik_e, prior_e, params_) + prop_e(params,params_) 
        denominator = post_e(lik_e, prior_e, params) + prop_e(params_,params)
        post_l.append((numerator, denominator)) 
        alpha = min(maxmin(numerator - denominator), 0)
        
        if np.isnan(alpha):
            print("!!!! ALPHA IS NAN !!!!")
            print(maxmin(numerator - denominator))
            print(lik_e(data, params_))
            print(params_)
            
        alphas.append(alpha)
        if np.random.uniform() < np.exp(alpha): 
            params = params_ 
            params_copy = copy.deepcopy(params)
            l.append(params_copy)
            idx.append(burnin+i)
            
    return l, idx, alphas, params_l, post_l

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
