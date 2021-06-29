"""
Match Pursuit Sparse Signal Recovery

References:
    Hammeed, Maxin Abdulrasool, "Comparative Analysis of Orthogonam Matching Pursuit
    and least angle regression," Michigan State University, A Thesis for Masters of
    Science, 2012.
"""

import numpy as np
import scipy.linalg as lin

# terminate with output signal has sparsity k
def sparsity_term(k, y, r, x_hat): 
    return int(np.linalg.norm(x_hat, 0, axis=0)) == k

# terminate when output signal has p percentage of signal
def percent_term(p, y, r, x_hat):
    y_l2p = np.linalg.norm(y) * (1.0-p)
    return y_l2p > np.linalg.norm(r)

# terminate when energy (L2 norm) of residual ie below residual
def epsilon_term(e, y, r, x_hat): 
    return e > np.linalg.norm(r)

def mp(y, A, term, param):
    """
    omp_sparsity: orthogonal match pursuit with configurable termination criteria

    Parameters:
        y: `compressed samples`
        A: `sampling matrix`
        term: `termination function from above`
        param: { k: `sparsity` | p: `percent` | e: `epsilon` }

    Returns:
        x_hat: `reconstructed signal`
    """

    r = y       # init residual
    x_hat = np.zeros((A.shape[1], 1)) # null output init

    while not term(param, y, r, x_hat):
        c = np.dot(A.T, r)               # correlation
        ind = np.argmax(np.abs(c))       # find abs support of max correlation
        x_hat[ind] += c[ind]             # update ind coefficient
        r = y - np.dot(A, x_hat)         # update residual

    return x_hat
