"""
Orthogonal Match Pursuit Sparse Signal Recovery

References:
    J. A. Tropp and A. C. Gilbert, "Signal Recovery From Random Measurements
    Via Orthogonal Matching Pursuit," in IEEE Transactions on Information Theory,
    vol. 53, no. 12, pp. 4655-4666, Dec. 2007.

    Hammeed, Maxin Abdulrasool, "Comparative Analysis of Orthogonal Matching Pursuit
    and least angle regression," Michigan State University, A Thesis for Masters of
    Science, 2012.
"""

import numpy as np
import scipy.linalg as lin

# terminate with output signal has sparsity k
def sparsity_term(k, y, r, x): 
    return int(np.linalg.norm(x, 0, axis=0)) == k

# terminate when output signal has p percentage of signal
def percent_term(p, y, r, x):
    y_l2p = np.linalg.norm(y) * (1.0-p)
    return y_l2p > np.linalg.norm(r)

# terminate when energy (L2 norm) of residual ie below residual
def epsilon_term(e, y, r, x): 
    return e > np.linalg.norm(r)

def omp(y, A, term, param):
    """
    orthogonal match pursuit with configurable termination criteria

    Parameters:
        y: `compressed samples`
        A: `sampling matrix`
        term: `termination function from above`
        param: { k: `sparsity` | p: `percent` | e: `epsilon` }

    Returns:
        x_hat: `reconstructed signal`
    """

    r = y       # init residual
    lamda = []  # list of support loc inds
    phi = []    # list of support vectors
    x = np.array([]) # null output init

    while not term(param, y, r, x):
        c = np.dot(A.T, r)               # correlation
        ind = np.argmax(np.abs(c))       # find abs support of max correlation

        lamda.append(ind)                # update support
        phi.append(A[:, ind:(ind+1)])    # update support locs 

        P = np.concatenate(phi, axis=1)  # compose the submatrix
        x = np.linalg.lstsq(P, y, rcond=None)[0]  # min ||y-P*x||_2
        r = y - np.dot(P, x)                      # update residual

    # embed coefficients in support
    n = A.shape[1]
    x_hat = np.zeros((n, 1))
    x_hat[lamda] = x

    return x_hat
