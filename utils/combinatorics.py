import numpy as np
# COMBINATORIAL FUNCTIONS 
def choose(n,k): 
    prod = 1
    for i in range(k): 
        prod *= (n-i)/(k-i)
    return prod

def log_choose(n,k):
    if n < k: return -700  
    if k==0 or k==n: return 0 
    prod = 0
    for i in range(k): 
        prod += (np.log(n-i) - np.log(k-i))
    return prod 

def factorial(n): 
    prod = 1 
    for i in range(n): 
        prod *= (n-i) 
    return prod 

def power(n,k): 
    prod = 1 
    for _ in range(k): 
        prod *= n
    return prod

# number of ways to fill n spaces with k distinct objects 
def fill_distinct(n,k): 
    s = 0 
    mult = 1 
    for i in range(k):
        s += mult * power(k-i, n) * choose(k,k-i) 
        mult *= -1
    return s 