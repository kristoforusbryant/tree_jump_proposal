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
    k = st.geom.rvs(geom_p)
    T_set = []
    for _ in range(k):  
        T_set.append(list(sample_tree(n).edges))
    E = set()
    for t in T_set: E = E.union(set(t))
    G = set_to_dic_list(n, E)
    return (G,T_set)

def prior_k_e(k, p=0.5): 
    return np.log(p) + (k-1) * np.log(1-p)

def prior_r_e(p0,p1, method=method): 
    maxmin = lambda x: max(min(x,np.exp(700)),np.exp(-700))
    d0,ts0 = p0 # previous
    d1,ts1 = p1 # proposed
    k0,k1 = (len(ts0), len(ts1))
    tau = maxmin(np.power(len(d0), (len(d0)-2)))
    
    if method == "MCMC": 
        mcmc0 = maxmin(np.log(Y(d0, k0, numIters=1000)))
        mcmc1 = maxmin(np.log(Y(d1, k1, numIters=1000)))
        approx = mcmc1 - mcmc0 
        prior_k_r = prior_k_e(k1, geom_p) - prior_k_e(k0, geom_p)

        ratio = approx + prior_k_r
    
    if method == "bounds": 
        lb = lower_bound(d1,d0)
        ub = upper_bound(d1,d0,ts1,ts0)

        approx = (ub + lb) / 2 # geometric average
        normalising = np.log(k0+1) - np.log(tau - k0) if k1 > k0 else np.log(tau - k0+1) - np.log(k0) 
        prior_k_r = prior_k_e(k1, geom_p) - prior_k_e(k0, geom_p)

        ratio = approx + normalising + prior_k_r

    return ratio 

# proposal  
def prop_s(params):
    G,T_set = copy.deepcopy(params)
    k = len(T_set)
    sample = list(sample_tree(n).edges)
    if sample in T_set and len(T_set) > 1: 
        T_set.remove(sample)
        E = set()
        for t in T_set: E = E.union(set(t))
        G = set_to_dic_list(n, E)
        for i in G.keys(): 
            G[i].sort() # sort in place 
        return (G,T_set)
    else: 
        T_set.append(sample)
        for i,j in sample:
            if i not in G[j]: 
                G[i].append(j)
                G[j].append(i)
        for i in G.keys(): 
            G[i].sort()
        return (G,T_set)

prop_e = lambda p: 0
prop_r_e = lambda p0,p1: prop_e(p1) - prop_e(p0)

# likelihood 
def lik_e(data, params): 
    maxmin = lambda x: max(min(x,np.exp(700)),np.exp(-700))
    U = lambda x: np.transpose(x) @ x
    D_star = constrained_cov(params, D + U(data), np.eye(D.shape[0]))
    log_prob_approx = laplace_approx(params, delta + len(data), D_star) - laplace_approx(params, delta, D) # as log prob 
    return log_prob_approx

lik_r_e = lambda data, p0, p1: lik_e(data, p1[0]) - lik_e(data, p0[0])