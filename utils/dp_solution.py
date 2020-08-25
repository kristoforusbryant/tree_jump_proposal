import numpy as np 
import networkx as nx
from .brute_force import count_spanning_partitions_brute

# Implement Dynamic Programming
def count_spanning_trees(G): 
    A = nx.adj_matrix(G).todense()
    D = np.diag(np.array(A.sum(axis=0))[0])
    L = D - A
    eigv = np.linalg.eigvals(L)
    return int(np.round(np.prod(eigv[eigv > 0.0000001])/L.shape[0]))

def count_spanning_partitions(G,k, memo={}):  
    #if in memo  
    if tuple(G.edges) in memo.keys(): 
        return memo[tuple(G.edges)], memo # is tuple the best hashable to be the key? consider 2^k encoding. 
    
    #if base cases
    if len(G.edges) == (len(G.nodes) + k - 2): # fix this conditions with the minimal cycle argument 
        res = count_spanning_partitions_brute(G,k)
        memo[tuple(G.edges)] = res
        return res, memo 
    
    #if nonsense 
    if len(G.edges) < len(G.nodes): 
        raise ValueError
    
    #else 
    tau = round(count_spanning_trees(G))
    total = 1
    for i in range(k): total *= (tau - i) / (k-i)
    
    sum_connected = 0 
    
    for i in range(1,np.power(2, len(G.edges)) -1):  # exclude the null and the graph G
        binary = (len(G.edges) - len(bin(i)[2:])) * '0' + bin(i)[2:]
        if sum([int(b) for b in binary]) < len(G.nodes): continue
        G_new = nx.Graph()
        G_new.add_nodes_from(list(G.nodes))
        G_new.add_edges_from([list(G.edges)[j] for j,b in enumerate(binary) if int(b)])
        if not nx.is_connected: print(G_new.edges)
        if nx.is_connected(G_new):
            G_next = G_new.copy()
            sums, memo = count_spanning_partitions(G_next, k, memo)
            sum_connected += sums
            
    memo[tuple(G.edges)] = total - sum_connected
    return round(total - sum_connected), memo