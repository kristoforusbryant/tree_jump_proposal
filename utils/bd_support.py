import numpy as np 

def complete_graph(n): 
    return {i:list(range(n)) for i in range(n)}
def complement_graph(d):
    d_ = complete_graph(len(d))
    for i in range(len(d)): 
        d_[i] = list(set(d_[i]) - set(d[i]) - set([i]))
    return d_
def one_connected_comp(d, start):
    to_visit = [start]
    visited = set()
    while to_visit:
        temp = set()
        for v in to_visit: 
            visited = visited.union({v})
            temp = temp.union(set(d[v]))
        to_visit = list(temp - visited)
    return visited 
def connected_comps(d): 
    s = list(range(len(d)))    
    acc = []
    while s: 
        cc = one_connected_comp(d,s[0])
        acc.append(cc)
        s = list(set(s) - cc)
    return acc 
def connection_graph(n1,n2,g):
    d = {}
    for v in list(n1):
        d[v] = list(set(g).intersection(n2))
    return d
## TODO: Check randomness
def random_edge(d): 
    k = np.random.choice(range(len(d)))
    return (k,np.random.choice(d[k]))

def reindex(d):
    idx = list(d.keys()) 
    new_idx = {idx[i]:i for i in range(len(d))}
    d_ ={new_idx[k]: [new_idx[i] for i in v] for k,v in d.items()}
    return d_

def subgraph(g, nodes):
    d = {}
    for n in list(nodes):
        d[n] = list(set(g[n]).intersection(set(nodes))) 
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
    
def lift_HG(hg, cc, g):
    lifted = []
    for k,l in hg.items(): 
        for v in l: 
            if k < v:
                lifted.append(connection_graph(cc[k], cc[v]), g)
    return lifted #graphs containing lifted edges 