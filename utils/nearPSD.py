import numpy as np
import copy 

def Projection_U(A, W):
    assert np.allclose(A, A.T, rtol=1e-5, atol=1e-8), "A not symmetric"
    assert np.allclose(W, W.T, rtol=1e-5, atol=1e-8), "W not symmetric"
    assert A.shape == W.shape, "A and W have different shapes"
    
    n = A.shape[0]
    try: 
        W_inv = np.linalg.inv(W)
    except: 
        print("W is singular")
        return -1
    theta = np.linalg.solve(W_inv * W_inv, np.diag(A) - np.ones(n))
    
    return A - (W_inv @ np.diag(theta) @ W_inv)

def Projection_S(A, W): 
    assert np.allclose(A, A.T, rtol=1e-5, atol=1e-8), "A not symmetric"
    assert np.allclose(W, W.T, rtol=1e-5, atol=1e-8), "W not symmetric"
    assert A.shape == W.shape, "A and W have different shapes"
    
    n = A.shape[0]
    try: 
        theta, Q = np.linalg.eig(W)
    except: 
        print("W is singular")
        return -1
    W_sqrt = Q @ np.diag(np.sqrt(theta)) @ np.transpose(Q)
    W_sqrt_inv = np.linalg.inv(W_sqrt)
    
    theta, Q = np.linalg.eig(W_sqrt @ A @ W_sqrt)
    D = np.diag(np.maximum(theta, np.zeros(n)))
    WAW_plus = Q @ D @ np.transpose(Q)

    return W_sqrt_inv @ WAW_plus @ W_sqrt_inv 

def weighted_norm(A, W): 
    try: 
        theta, Q = np.linalg.eig(W)
    except: 
        print("W is singular")
        return -1
    W_sqrt = Q @ np.diag(np.sqrt(theta)) @ np.transpose(Q)
    A_ = W_sqrt @ A @ W_sqrt
    
    return np.sum(A_ * A_)

def nearPSD(A, W): 
    # W is positive definite matrix having the same shape as A 
    A_ = np.zeros(A.shape) - 10000
    dS = np.zeros(A.shape)
    
    count = 0 
    while (weighted_norm(A - A_,W) > 1e-10):
        count +=1 
        A_ = A
        Rk = A - dS 
        Xk = Projection_S(Rk, W)
        dS = Xk - Rk
        A = Projection_U(Xk, W)
        
    # one last projection to S 
    Rk = A - dS
    A = Projection_S(Rk, W)

    return (A)