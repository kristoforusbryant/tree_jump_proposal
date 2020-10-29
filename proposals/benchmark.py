### BENCHMARK 
# prior
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
prop_e = lambda p: 0
prop_r_e = lambda p0,p1: prop_e(p1) - prop_e(p0)

# likelihood 
def lik_e(data, params): 
    U = lambda x: (np.transpose(x) @ x).round(2)
    D_star = constrained_cov(params, D + U(data), np.eye(D.shape[0]))
    log_prob_approx = laplace_approx(params, delta + len(data), D_star) - laplace_approx(params, delta, D) # as log prob 
    return log_prob_approx

lik_r_e = lambda data, p0, p1: lik_e(data, p1) - lik_e(data, p0)