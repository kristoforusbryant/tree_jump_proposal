import numpy as np 
import networkx as nx 

# Function from Prufer Sequences to Trees
def prufer_to_trees(a): 
    # array of n integers from 0 to n+1 (inclusive)
    n = len(a)
    T = nx.Graph()
    T.add_nodes_from(list(range(n+2)))
    deg = [1 for _ in range(len(T.nodes))]
    for i in a: 
        deg[i] += 1
        
    for i in a: 
        for j in range(len(T.nodes)): 
            if deg[j] == 1: 
                T.add_edge(i,j)
                deg[j] -= 1
                deg[i] -= 1
                break
    
    u, v = 0, 0
    for i in range(len(T.nodes)): 
        if deg[i] == 1:
            if u == 0: 
                u = i
            else: 
                v = i 
                break 
    
    T.add_edge(u,v)
    deg[u] -= 1
    deg[v] -= 1
    
    return T

def sample_tree(n): 
    seq = [np.random.randint(0,n) for _ in range(n-2)]
    return(prufer_to_trees(np.array(seq)))