#!/usr/bin/python

"""
Match Pursuit Sparse Signal Recovery

References:
    Hammeed, Maxin Abdulrasool, "Comparative Analysis of Orthogonal Matching Pursuit
    and least angle regression," Michigan State University, A Thesis for Masters of
    Science, 2012.
"""

import numpy as np
import scipy.linalg as lin

def mp(y, A, term, param):
    """
    omp_sparsity: orthogonal match pursuit with configurable termination criteria

    Parameters:
        y: `compressed samples`
        A: `sampling matrix`
        term: `termination function`
        param: `termination parameter` { k: `sparsity` | p: `percent` | e: `epsilon` }

    Returns:
        x_hat: `reconstructed signal`
    """

    r = np.copy(y)                       # init residual
    x_hat = np.zeros((A.shape[1], 1))    # null output init

    while not term(param, y, r, x_hat):
        c = np.dot(A.T, r)               # correlation
        ind = np.argmax(np.abs(c))       # find abs support of max correlation
        x_hat[ind] += c[ind]             # update ind coefficient
        r = y - np.dot(A, x_hat)         # update residual

    return x_hat

if __name__ == "__main__":
    import harness
    harness.test(mp)
    harness.test(mp, True)
