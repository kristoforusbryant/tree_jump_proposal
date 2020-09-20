import numpy as np 
import networkx as nx 
from itertools import combinations, chain
from .combinatorics import *

# Defining a Graph 
class Graph(): 
    def __init__(self, n, E):
        # Assume that n is a positive integer 
        # Assume that E is a set of tuples containing integers (increasing order) between 0 and n-1 (inclusive)  
        # TODO: write checks for these 
        self.n = n
        self.E = [(min(e), max(e)) for e in E] 
        self.IN = set()
        self.OUT = set()

    def draw(self):     
        G = nx.Graph()
        G.add_edges_from(self.E)
        nx.draw_shell(G, with_labels=True)
        plt.show()
    
    def edge_dict(self): 
        dic = { i:[] for i in range(self.n)} 
        for k,v in self.E: 
            dic[k].append(v)
            dic[v].append(k)
        return dic

def find_fundamental_cycle(T, e):
    e = (min(e), max(e))
    # maintain a queue of paths
    dic = T.edge_dict()
    queue = []
    # push the first path into the queue
    queue.append([e[1]])
    while queue:
        # get the first path from the queue
        path = queue.pop(0)
        # get the last node from the path
        node = path[-1]
        # path found
        if node == e[0]:
            path += [path[0]]
            path = [(min(path[i:i+2]), max(path[i:i+2])) for i in range(len(path)-1)]
            return path  # returns edges outlining the cycle
        # enumerate all adjacent nodes, construct a new path and push it into the queue
        for adjacent in dic[node]:
            new_path = list(path)
            new_path.append(adjacent)
            queue.append(new_path)

# Implementation of (Kapoor et al., 1992)
def all_spanning_trees(G): 
    # Preprocess 
    init_T = Graph(len(G.nodes), list(nx.minimum_spanning_tree(G).edges))
    init_T.E = init_T.E
    E_G = [(min(e), max(e)) for e in list(G.edges)]
    acc = [init_T.E]
    temp1 = [init_T]
    temp2 = []
    
    while temp1: 
        for T in temp1: 
            E_Tcomp = set(E_G) - set(T.E) - T.OUT 
            if not E_Tcomp: continue
            added = []
            for to_add in E_Tcomp:
                fc = set(find_fundamental_cycle(T, to_add)) - T.IN - set([to_add])
                if not fc: continue 
                removed = []
                for to_remove in fc:                    
                    E_new = T.E + [to_add]
                    if to_remove not in E_new: to_remove = (to_remove[1], to_remove[0]) 
                    E_new.remove(to_remove)
                    T_new = Graph(T.n, E_new)
                    T_new.IN = T.IN.union(set([to_add] + removed))
                    T_new.OUT = T.OUT.union(set([to_remove] + added))
                    acc.append(T_new.E)
                    temp2.append(T_new)
                    removed.append(to_remove)
                added.append(to_add)
        temp1 = temp2 
        temp2 = [] 
        
    return acc 

# Implement Brute Force Method
def count_spanning_partitions_brute(G, k):
    all_trees = all_spanning_trees(G)
    truth_count = 0
    set_E = set([(min(e), max(e)) for e in G.edges])
    
    for j in range(1,k+1): 
        for idxs in list(combinations(range(len(all_trees)), j)):
            truth_count += (set(chain.from_iterable([all_trees[i] for i in idxs])) == set_E) 
    return truth_count
