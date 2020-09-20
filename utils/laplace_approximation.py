import numpy as np
import copy
import random  

def BronKerbosch2(G, P, R=None, X=None):
    """ 
    Bron-Kerbosch Algorithm with Pivot from Bron and Kerbosch(1973)
    G: Graph as dict of lists
    """
    P = set(P)
    R = set() if R is None else R
    X = set() if X is None else X
    if not P and not X:
        yield R
    try:
        u = list(P.union(X))[0]
        S = P.difference(G[u])
    # if union of P and X is empty
    except IndexError:
        S = P
    for v in S:
        yield from BronKerbosch2(G, 
            P=P.intersection(G[v]), R=R.union([v]), X=X.intersection(G[v]))
        P.remove(v)
        X.add(v)
    
def mode(G, delta, D, N=50): 
    """
    G: Graph as dict of lists 
    """
    L = D / (delta - 2)
    K = np.eye(D.shape[0])
    C_l = list(BronKerbosch2(G, G.keys()))
    for i in range(N):
        K_ = copy.deepcopy(K)
        for c in C_l:
            not_c = tuple(sorted(list(set(range(len(G))) - set(c))))
            c = tuple(c)
            B = K[not_c,:][:, c]
            C = K[c,:][:, not_c]
            D = K[not_c,:][:, not_c]
            K[np.ix_(c, c)]= np.linalg.inv(L[c,:][:,c]) + (C @ np.linalg.inv(D) @ B) 
        if np.max(K - K_) < 1e-8: break 
    return K 

def constrained_cov(G, L, M, N=50):
    """
    The second cyclic algorithm in (Speed and Kiiveri, 1986)
    """
    C_l = list(BronKerbosch2(G, G.keys()))
    K = np.linalg.inv(M)
    K_ = K
    for i in range(N):
        for c in C_l:
            c = tuple(c)
            Q_inv = np.linalg.inv(K_) 
            Q_inv[np.ix_(c, c)] += np.linalg.inv(L[c,:][:,c]) - np.linalg.inv(K_[c,:][:,c]) # subset first, then take inv     
            K_ = np.linalg.inv(Q_inv)
        if np.max(K - K_) < 1e-8: break
        K = K_
    return K 

def laplace_approx(G, delta, D): 
    """
    Laplace Approximation as outlined by (Lenkoski and Dobra, 2011)
    """
    K = mode(G, delta, D)
    V = []
    for k,l in G.items():
        V.append((k,k))
        for v in l: 
            if k < v: V.append((k,v))
    h = -.5 * (np.trace(np.transpose(K) @ D) - (delta - 2) * np.log(np.linalg.det(K))) 
    H = hessian(K, V, delta, D)
    return np.exp(h) * np.power(2*np.pi, len(V)/2) * np.power(np.linalg.det(-H), -1/2)
    
def hessian(K, V, delta, D): 
    H = np.zeros((len(V), len(V)))
    K_inv = np.linalg.inv(K)
    for a in range(len(V)):
        i,j = V[a]
        one_ij = np.zeros(K_inv.shape)
        one_ij[(i,j),(j,i)] = 1 
        for b in range(len(V)): 
            k,l = V[b]
            one_kl = np.zeros(K_inv.shape)
            one_kl[(k,l),(l,k)] = 1
            H[a,b] = -.5*(delta-2)*np.trace(K_inv @ one_ij @ K_inv @ one_kl)
    return H