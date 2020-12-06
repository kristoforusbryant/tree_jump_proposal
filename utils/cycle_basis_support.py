import copy  
from .bd_support import subgraph 

def find_cycle(d, now, prev, visited):
    for v in d[now]: 
        if v == prev: continue
        if v in visited:
            return visited 
        else: 
            cycle = find_cycle(d, v, now, visited + [v])
            if cycle: 
                return cycle  
    else: 
        return None 
    
def cycle_basis(T): 
    # TODO: check if T is a spanning tree of its nodes
    basis = []
    for i in range(len(T)): 
        for j in range(i):
            if j in T[i]: continue
            T_ = copy.deepcopy(T)
            T_[i].append(j)
            T_[j].append(i)
            basis.append(subgraph(T_, find_cycle(T_, i, -1, [])))
    return basis

def basis_addition(b0,b1):  
    # TODO: consider whether using sets is faster
    g = copy.deepcopy(b0)
    for k,l in b1.items(): 
        for v in l: 
            if v < k: continue
            if v in g[k]: 
                g[k].remove(v)
                g[v].remove(k)
            else: 
                g[k].append(v)
                g[v].append(k)
    return g