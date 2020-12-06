import scipy.stats as st
from utils.laplace_approximation import laplace_approx, constrained_cov
from utils.sampling import sample_tree
from utils.prior_MCMC import Y
from utils.bounds import lower_bound, upper_bound

print("geom_p: " + str(geom_p)) 
print("method: " + str(method))

## TREE-PROPOSAL
def set_to_dic_list(n, E): 
    dic = {i:[] for i in range(n)}
    for i,j in E: 
        dic[i].append(j)
        dic[j].append(i)
    for i in dic.keys(): 
        dic[i].sort()
    return dic 

# prior
def prior_s():
    k = 1 
    T_set = []
    for _ in range(k):  
        T_set.append(list(sample_tree(n).edges))
    E = set()
    for t in T_set: E = E.union(set(t))
    G = set_to_dic_list(n, E)
    return (G,T_set)

def prior_r_e(p0,p1, method=method): 
    maxmin = lambda x: max(min(x,np.exp(700)),np.exp(-700))
    d0,ts0 = copy.deepcopy(p0) # previous
    d1,ts1 = copy.deepcopy(p1) # proposed
    k0,k1 = (1,1)

    mcmc0 = maxmin(np.log(Y(d0, k0, numIters=5000)))
    mcmc1 = maxmin(np.log(Y(d1, k1, numIters=5000)))
    approx = mcmc1 - mcmc0 

    ratio = approx 

    return ratio 

# proposal  
def prop_s(params):
    G,T_set = copy.deepcopy(params)
    k = len(T_set)

    sample = list(sample_tree(n).edges)
    while sample in T_set: 
        sample = list(sample_tree(n).edges)

    T_set = []
    G_new = {k:[] for k in range(n)}
    T_set.append(sample)
    for i,j in sample:
        if i not in G_new[j]: 
            G_new[i].append(j)
            G_new[j].append(i)
    for i in G_new.keys(): 
        G_new[i].sort()
    return (G_new,T_set)

prop_e = lambda p: 0
prop_r_e = lambda p0,p1: prop_e(p1) - prop_e(p0)

# likelihood 
def lik_e(data, params): 
    params = copy.deepcopy(params)
    maxmin = lambda x: max(min(x,np.exp(700)),np.exp(-700))
    U = lambda x: np.transpose(x) @ x
    D_star = constrained_cov(params, D + U(data), np.eye(D.shape[0]))
    log_prob_approx = laplace_approx(params, delta + len(data), D_star) - laplace_approx(params, delta, D) # as log prob 
    return log_prob_approx

lik_r_e = lambda data, p0, p1: lik_e(data, p1[0]) - lik_e(data, p0[0])