import scipy.stats as st
from utils.laplace_approximation import laplace_approx, constrained_cov
from utils.sampling import sample_tree
from utils.prior_MCMC import Y, sampleTreeFromConnected
from utils.bounds import lower_bound, upper_bound
from utils.bd_support import *
from utils.dp_solution import count_spanning_trees


geom_p = .5 
method = 'bounds'
birth_p = .5 
beta = .8

print("geom_p: " + str(geom_p)) 
print("method: " + str(method))
print("birth_p: " + str(birth_p)) 
print("beta: " + str(beta)) # probability of being smart   

def count_spanning_trees(G): 
    G = nx.from_dict_of_lists(G)
    A = nx.adj_matrix(G).todense()
    D = np.diag(np.array(A.sum(axis=0))[0])
    L = D - A
    eigv = np.linalg.eigvals(L)
    return int(np.round(np.prod(eigv[np.abs(eigv) > 0.0000001])/L.shape[0]))

# Smart Birth Death Proposal with Mixing 
## According to (Abhinav, Week 10) 

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

def prior_r_e(p0,p1, method = method): 
    maxmin = lambda x: max(min(x,np.exp(700)),np.exp(-700))
    d0,ts0 = p0 # previous
    d1,ts1 = p1 # proposed
    k0,k1 = (len(ts0), len(ts1))
    
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

def _prior_r_e(p0,p1): 
    maxmin = lambda x: max(min(x,np.exp(700)),np.exp(-700))
    d0,ts0 = p0 # previous
    d1,ts1 = p1 # proposed
    k0,k1 = (len(ts0), len(ts1))
    tau = maxmin(np.power(len(d0), (len(d0)-2)))

    lb = lower_bound(d1,d0)
    ub = upper_bound(d1,d0,ts1,ts0)

    approx = (ub + lb) / 2 # geometric average
    normalising = np.log(k0+1) - np.log(tau - k0) if k1 > k0 else np.log(tau - k0+1) - np.log(k0) 
    prior_k_r = prior_k_e(k1, geom_p) - prior_k_e(k0, geom_p)

    ratio = approx + normalising + prior_k_r

    return ratio 

# BD proposal 
def death(G,T_set): 
    sample = T_set[np.random.choice(len(T_set))]
    T_set.remove(sample)
    E = set()
    for t in T_set: E = E.union(set(t))
    G = set_to_dic_list(n, E)
    for i in G.keys(): 
        G[i].sort()
    return ((G,T_set), -np.log(len(T_set))) 
death_prior = lambda G,T_set: np.log(len(T_set))

def birth(G,T_set): 
    sample = T_set[0]
    while sample in T_set: sample = list(sample_tree(n).edges)
    T_set.append(sample)
    for i,j in sample:
        if i not in G[j]: 
            G[i].append(j)
            G[j].append(i)
    for i in G.keys(): 
        G[i].sort()
    return ((G,T_set), count_spanning_trees(G) - len(T_set)) 
birth_prior = lambda G,T_set: count_spanning_trees(G) - len(T_set)

# Smart-Birth BD proposal 
def gen_tree_smart(G): 
    G_ = complement_graph(G) 
    C = connected_comps(G_)
    Tc = []
    lifted = []
    lp = 0 
    for c in C: 
        S = subgraph(G,c)
        Tc += as_edge_list(sampleTreeFromConnected(S))
        lp -= np.log(count_spanning_trees(S))
    for i in range(len(C)): 
        for j in range(i):
            S = lift_HG(C[i], C[j], G)
            lifted.append(random_edge(S))
            lp -= np.log(n_edges(S)) 
    return Tc + lifted, lp 
    
def smart_birth(G,T_set, beta):
    if np.random.uniform() < beta: 
        sample, lp = gen_tree_smart(G)
        T_set.append(sample)
        for i,j in sample:
            if i not in G[j]: 
                G[i].append(j)
                G[j].append(i)
        for i in G.keys(): 
            G[i].sort()
        return ((G,T_set), lp + np.log(beta))
    else: 
        params, lp = birth(G, T_set)
        return (params, lp + np.log(1-beta))

def smart_birth_p(G,T_set, beta): # computing probability needs requires one pass through the algorithm  
    if np.random.uniform() < beta: 
        sample, lp = gen_tree_smart(G)
        T_set.append(sample)
        for i,j in sample:
            if i not in G[j]: 
                G[i].append(j)
                G[j].append(i)
        for i in G.keys(): 
            G[i].sort()
        return lp + np.log(beta)
    else: 
        params, lp = birth(G, T_set)
        return lp + np.log(1-beta)
    
def bd_prop(params, p=birth_p, beta=beta):
    G,T_set = copy.deepcopy(params)
    if len(T_set) < 2: p = 1
    lp = 0
    # Birth or Death 
    if np.random.uniform() < p or len(T_set) == 1:
        params, lp_p1_p0 = smart_birth(G,T_set, beta) # logP(p0 | p1) 
        lp_p0_p1 = death_p(*params) 
        lp = lp + np.log(1-p) + lp_p0_p1
        lp = lp - np.log(p) - lp_p1_p0
    else: 
        params, lp_p1_p0 = death(G,T_set)
        lp_p0_p1 = smart_birth_p(G,T_set, beta)
        lp = lp + np.log(p) + lp_p0_p1
        lp = lp - np.log(1-p) - lp_p1_p0
    # Mixing 
    ## small enough to ignore
    return (death(*(birth(*params)[0]))[0], lp)

def prop_e(p1,p0,p=birth_p): # q(p1 | p0) 
    n = len(p0[0])
    T0 = p0[1]
    T1 = p1[1] 
    assert len(T0) == (len(T1) - 1) or len(T0) == (len(T1) + 1), "Graphs doesn't differ by exactly one tree" 
    
    if len(T1) > len(T0): # birth event
        return np.log(p) -  np.log(np.power(n,n-2) - len(T0))
    else: # death event 
        return np.log(1-p) - np.log(len(T0))
    
def prop_r_e(p0,p1): return prop_e(p0, p1) - prop_e(p1, p0)    

prop = bd_prop