from utils.combinatorics import choose

### BENCHMARK 
# prior
if PRIOR_TYPE=='uniform': 
    def prior_s():
        G = {i:[] for i in range(n)}
        for i in range(n): 
            for j in range(i):
                if np.random.uniform() > 0.5:
                    G[i].append(j)
                    G[j].append(i)
        for i in G.keys(): 
            G[i].sort()
        return G
    prior_e = lambda p: 0
    prior_r_e = lambda p0,p1: prior_e(p1) - prior_e(p0) 
    
elif PRIOR_TYPE=='size-based': 
    def prior_s():
        G = {i:[] for i in range(n)}
        for i in range(n): 
            for j in range(i):
                if np.random.uniform() > 0.5:
                    G[i].append(j)
                    G[j].append(i)
        for i in G.keys(): 
            G[i].sort()
        return G
    def prior_e(p): 
        m = int(n*(n-1)/2)
        k = n_edges(p)
        return( - np.log(m+1) - np.log(choose(m, k)) ) 
    prior_r_e = lambda p0,p1: prior_e(p1) - prior_e(p0) 

def prop_s(G): 
    i,j = random.choice(list(combinations(list(G.keys()), 2)))
    G = copy.deepcopy(G)
    if j in G[i]: 
        G[i].remove(j)
        G[j].remove(i)
    else: 
        G[i].append(j)
        G[i].sort() # is sorting necessary?
        G[j].append(i) 
        G[j].sort()
    return G 
prop_e = lambda p0,p1: 0
prop_r_e = lambda p0,p1: prop_e(p1) - prop_e(p0)

# likelihood 
def lik_e(data, params): 
    U = lambda x: (np.transpose(x) @ x).round(2)
    D_star = constrained_cov(params, D + U(data), np.eye(D.shape[0]))
    log_prob_approx = laplace_approx(params, delta + len(data), D_star) - laplace_approx(params, delta, D) # as log prob 
    return log_prob_approx

def lik_r_e(data, p0, p1): 
    return lik_e(data, p1) - lik_e(data, p0)