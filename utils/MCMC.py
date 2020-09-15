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
    
    post_e = lambda lik,prior,params: lik(data, params) * np.prod(prior(params))
    # Initialisation 
    while True: 
        params = prior_s()
        if post_e(lik_e, prior_e, params) != 0: 
            break
    # Burn in 
    for _ in tqdm(range(burnin)):
        params_ = prop_s(params)
        numerator = (post_e(lik_e, prior_e, params_) * prop_e(params,params_)) 
        denominator = (post_e(lik_e, prior_e, params) * prop_e(params_,params))
        
        if denominator > 0: 
            alpha = min(numerator / denominator, 1)
        else: 
            alpha = 1
        if np.random.uniform() < alpha: 
            params = params_ 
            params_copy = copy.deepcopy(params)
    # Sampler
    for _ in tqdm(range(N)):
        params_l.append(params)
        params_ = prop_s(params)
        numerator =  (post_e(lik_e, prior_e, params_) * prop_e(params,params_))
        denominator = (post_e(lik_e, prior_e, params) * prop_e(params_,params))
        if denominator > 0: 
            alpha = min(numerator / denominator, 1)
        else: 
            alpha = 1
        alphas.append(alpha)
        if np.random.uniform() < alpha: 
            params = params_ 
            params_copy = copy.deepcopy(params)
            l.append(params_copy)
    return l, alphas, params_l

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
        sample = []
        for i in range(len(self.param_l)): 
            sample.append(self.fun_l[i](data[i], *self.param_l[i]))
        return tuple(sample)
