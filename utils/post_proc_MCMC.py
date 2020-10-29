import numpy as np 
import matplotlib.pyplot as plt 

def create_edge_matrix(samples): 
    n = len(samples[0])
    edge_M = np.zeros((n,n))
    for d in samples:
        for k, li in d.items(): 
            for v in li: 
                if k > v:
                    edge_M[k,v] +=1
                    edge_M[v,k] +=1
    edge_M = edge_M/len(samples)
    return edge_M

def heatmap(matrix, save=None): 
    plt.figure(figsize=(10, 10))
    plt.imshow(matrix, cmap='viridis')
    plt.colorbar()    
    if save: 
        plt.savefig(save)
    plt.show()