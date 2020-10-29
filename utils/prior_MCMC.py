import numpy as np 

def sampleTreeFromConnected(d): 
    """
    Aldous-Broder Algorithm
    sampling scheme for uniform tree from connected graph 
    """
    n = len(d)
    visited = []
    x = np.random.choice(list(d.keys()))
    T = {k:[] for k in d}
    
    # init
    visited.append(x)
    
    # random walk on the graph  
    while len(visited) < n: 
        nbh = d[x]
        r = np.random.choice(nbh)
        if (r not in visited): 
            visited.append(r) 
            T[r].append(x)
            T[x].append(r) 
            x = r 
        else: 
            x = r
    return T

def k_core(d, K):
    idx_to_remove = set()
    # do once 
    for k, l in d.items():
        if len(l) < K: 
            idx_to_remove = idx_to_remove.union({k}) 
    for k in idx_to_remove: del d[k]
    
    while idx_to_remove: 
        # remove keys
        for k in d:
            for i in list(idx_to_remove):
                if i in d[k]: d[k].remove(i)
        idx_to_remove = set() # reset 
        
        for k, l in d.items(): 
            if len(l) < K: 
                idx_to_remove.union({k})
        for k in idx_to_remove: del d[k]
    return d

def n_edges(d): 
    count = 0 
    for k,l in d.items(): 
        for v in l: 
            if k < v: 
                count += 1
    return count

def add_graph(d0,d1): 
    for k, l in d1.items(): 
        d0[k] = list(set(d0[k]).union(set(d1[k])))
    return d0

def _Y(d, k, numIters):
    success = 0
    for it in range(numIters): 
        T = {k:[] for k in d}
        for i in range(k):
            T = add_graph(T, sampleTreeFromConnected(d))
        for k in T: T[k].sort()
        if T == d: 
            success += 1
    return success

def Y(d, k, numIters=1000, numWorkers=1):
    # basic bounds checking
    n, m = len(d), n_edges(d)
    if (k*(n-1) < m):
        return 0.0
    
    # reduce to 2-core
    G2core = k_core(d,2)
    n2core = len(G2core)
    if (n2core == 0):
        return 1
    
    if (numWorkers > 1):
        batch = numIters // numWorkers
        pool = mp.Pool(numWorkers)
        result = [pool.apply(_Y, args=(G2core, k, batch)) for i in range(numWorkers)]
        pool.close()
        success = sum(result)
        return success/(numWorkers*batch)
    else:
        return _Y(G2core, k, numIters)/numIters