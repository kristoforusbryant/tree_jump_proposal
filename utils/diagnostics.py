import numpy as np 

def thin(a, lag=1):
    indices = [i for i in range(len(a)) if i % lag == 0]
    return a[indices]

def AC(a, lag=1):
    M = len(a)
    a = a - np.mean(a)
    c_f_0 = np.sum(a[:M] * a[:M]) / M
    return np.sum(a[:M - lag] * a[lag:M]) / (M-lag) / c_f_0

def ACF(a, lag=1):
    return [AC(a[:i], lag=lag) for i in range(1,len(a))]

def AIC_time(a, M = None):
    # Implementation according to https://emcee.readthedocs.io/en/stable/tutorials/autocorr/
    if not M: M = len(a)
    a = a - np.mean(a)
    c_f_0 = np.sum(a[:M] * a[:M]) / M
    return 1 + 2* np.sum([ np.sum(a[:M - i] * a[i:M]) / (M-i) / c_f_0 for i in range(M)])