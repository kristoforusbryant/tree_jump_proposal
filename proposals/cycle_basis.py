from utils.cycle_basis_support import * 
from utils.combinatorics import choose

### CYCLE-BASIS
# prior (given the basis)
if PRIOR_TYPE=='uniform': 
    def prior_s(basis=basis): 
        G = {i:[] for i in range(n)}
        b_index = []
        for i in range(len(basis)): 
            if np.random.uniform() > 0.5:
                G = basis_addition(G, basis[i])
                b_index.append(i)
        return (G, b_index) 
    prior_e = lambda p: 0
    prior_r_e = lambda p0,p1: prior_e(p1) - prior_e(p0) 
    
elif PRIOR_TYPE=='size-based': 
    def prior_s(basis=basis): 
        G = {i:[] for i in range(n)}
        b_index = []
        for i in range(len(basis)): 
            if np.random.uniform() > 0.5:
                G = basis_addition(G, basis[i])
                b_index.append(i)
        return (G, b_index) 
    def prior_e(p): 
        G, b_ind = p 
        m = len(basis) 
        k = len(b_ind)
        return(- np.log(m+1) - np.log(choose(m, k)))
    prior_r_e = lambda p0,p1: prior_e(p1) - prior_e(p0) 
    
def prop_s(params):
    G, b_i = params
    i = np.random.choice(len(basis))
    b_i_ = copy.deepcopy(b_i)
    if i in b_i:
        b_i_.remove(i)
    else: 
        b_i_.append(i)
    return (basis_addition(G, basis[i]), b_i_)
prop_e = lambda p0, p1: 0
prop_r_e = lambda p0,p1: prop_e(p0,p1) - prop_e(p1, p0)

# likelihood 
def lik_e(data, params): 
    params = params[0]
    maxmin = lambda x: max(min(x,np.exp(700)),np.exp(-700))
    U = lambda x: np.transpose(x) @ x
    D_star = constrained_cov(params, D + U(data), np.eye(D.shape[0]))
    log_prob_approx = laplace_approx(params, delta + len(data), D_star) - laplace_approx(params, delta, D) # as log prob 
    return log_prob_approx

def lik_r_e(data, p0, p1): 
    return lik_e(data, p1) - lik_e(data, p0)