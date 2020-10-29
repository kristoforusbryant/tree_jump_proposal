import numpy as np 
from utils.combinatorics import choose, log_choose
import copy

def set_to_dict_lists(s, n): 
    dic = {i:[] for i in range(n)}
    for e in s: 
        if e[1] not in dic[e[0]]: 
            dic[e[0]].append(e[1])
            dic[e[1]].append(e[0]) 
    return dic

def compliment(d): 
    dic = {}
    s = set(d.keys())
    for k, vals in d.items(): 
        dic[k] = list(s - set(vals) - set([k]))
    return dic

def add_graph(d1,d2): 
    d = copy.deepcopy(d1)
    for k,v in d2.items():
        d[k] = list(set(d1[k]).union(set(v)))
    return d

def w_dict_to_adj(d): 
    A = np.zeros((len(d), len(d)))
    for k,l in d.items(): 
        for v in l: 
            A[k, v[0]] = v[1]
    return A 

def dict_to_adj(d): 
    A = np.zeros((len(d), len(d)))
    for k,l in d.items(): 
        for v in l: 
            A[k, v] = 1
    return A 

def difference(d1,d2): 
    dic = {}
    for i in d1.keys(): 
        dic[i] = list(set(d2[i]) - set(d1[i])) 
    return dic

def dfs(d, temp, v, visited):
    visited[v] = True 
    temp.append(v)
    for i in d[v]:
        if not visited[i]: 
            temp = dfs(d, temp, i, visited)
    return temp

def find_connected(d):
    visited = {i:False for i in d.keys()}
    to_visit = set(d.keys())
    components = []
    while to_visit:
        c = dfs(d, [], to_visit.pop(), visited)
        to_visit = to_visit - set(c)
        components.append(c)
    return components 

def weighted_quotient_graph(d1, d2):
    ddiff = difference(d1,d2)
    components = find_connected(ddiff)
    n = len(components)
    q = {i:[] for i in range(len(components))}
    for i in range(len(components)):
        for j in range(i): 
            w = 0
            for k in components[i]:
                w += len(set(d2[k]).intersection(set(components[j]))) 
            if w > 0: 
                q[i].append((j,w))
                q[j].append((i,w))
    return q

def count_bridging_trees(d1,d2): 
    q = weighted_quotient_graph(d1,d2)
    A = w_dict_to_adj(q)
    L = np.diag(A.sum(1)) - A 
    return round(np.linalg.det(L[1:,1:])) # any minor of L 

def count_spanning_trees(d):
    A = dict_to_adj(d)
    L = np.diag(A.sum(1)) - A 
    return round(np.linalg.det(L[1:,1:])) # any minor of L 

def lower_bound(d1,d0, as_log=True): 
    "Ratio of G1/G0"
    m1 = sum([len(d1[i]) for i in range(len(d1))])
    m0 = sum([len(d0[i]) for i in range(len(d0))])
    
    if m0 > m1: 
        if as_log:
            return -np.log(count_bridging_trees(d1,d0))
        else:
            return 1/count_bridging_trees(d1,d0)
    elif m1 > m0: 
        if as_log: 
            return np.log(count_bridging_trees(d0,d1))
        else: 
            return count_bridging_trees(d0,d1)
    #TODO: handle the case where graph doesn't change from addition/removal of trees 
    else: return -700 # if graph doesn't change => do not accept proposal 

# transform this to log scale 
def tree_count_product(tree_seq, as_log=True):
    G = tree_seq[0]
    G_new = tree_seq[0]
    prod = 0 if as_log else 1 
    for i in range(1, len(tree_seq)): 
        for k,v in tree_seq[i].items(): 
            G_new[k] = list(set(G[k]).union(set(v)))
        if as_log: 
            prod += np.log(count_bridging_trees(G, G_new))
        else: 
            prod *= count_bridging_trees(G, G_new)
        G = copy.deepcopy(G_new)
    return prod 

def upper_bound(d1,d0,ts1,ts0, as_log=True): 
    "Ratio of G1/G0"
    assert np.abs(len(ts0) - len(ts1)) == 1, 'proposal k does not differ by one'
    assert len(d1) == len(d0), 'number of nodes does not match'    
    
    t1_dic_list = [set_to_dict_lists(t, len(d1)) for t in ts1] 
    if len(ts0) > len(ts1): 
        tau = count_spanning_trees(d0)
        if as_log: 
            return tree_count_product(t1_dic_list) - log_choose(tau,len(t1_dic_list)+1)
        else: 
            return tree_count_product(t1_dic_list, as_log=False)/choose(tau,len(t1_dic_list)+1)
    else: 
        tau = count_spanning_trees(d1)
        if as_log: 
            return log_choose(tau,len(t1_dic_list)) - tree_count_product(t1_dic_list[:-1])
        else: 
            return choose(tau,len(t1_dic_list))/tree_count_product(t1_dic_list[:-1], as_log=False)

# def upper_bound(d1,d2,tree_seq): 
#     "Ratio of G1/G2, where tree_seq is the trees that builds G1"
#     m1 = sum([len(d1[i]) for i in range(len(d1))])
#     m2 = sum([len(d2[i]) for i in range(len(d2))])
    
#     if m2 > m1: 
#         tau = count_spanning_trees(d2)
#         return tree_count_product(tree_seq)/choose(tau,len(tree_seq)+1)
#     elif m1 > m2: 
#         tau = count_spanning_trees(d1)
#         return choose(tau,len(tree_seq))/tree_count_product(tree_seq[:-1])
#     else: return -1