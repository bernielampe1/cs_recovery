"""
Stagewise Orthogonal Match Pursuit Sparse Signal Recovery

References:
    D. Donoho, Y. Tsaig, I. Drori, J Starck, "Sparse Solution of Underdetermined Linear Equations
    by Stagewise Orthogonal Matching Pursuit,", March 2006

    Hammeed, Maxin Abdulrasool, "Comparative Analysis of Orthogonal Matching Pursuit
    and least angle regression," Michigan State University, A Thesis for Masters of
    Science, 2012.

    Angshul Majumdar, "Compressed Sensing for Engineers", CRC Press, 3 Dec 2018.
"""

import numpy as np
import scipy.linalg as lin

def stomp(y, A, sigma, N):
    """
    stagewise orthogonal match pursuit with configurable termination criteria

    Parameters:
        y: `compressed samples`
        A: `sampling matrix`
        sigma: `number of std deviations (2-3)`
        N: `number of iterations`

    Returns:
        x_hat: `reconstructed signal`
    """

    r = np.copy(y)       # init residual
    lamda = set()        # set of support indicies
    x = np.array([])     # null output init
    n = A.shape[1]       # size of solution

    for t in range(N):
        c = np.dot(A.T, r)                  # correlation
        thr = sigma * np.std(r)             # compute number of std devs for threshold
        inds = np.nonzero(np.abs(c) > thr)  # find abs support above threshold
        lamda.update(inds[0].tolist())      # update support

        P = A[:, list(lamda)]                     # compose the submatrix
        x = np.linalg.lstsq(P, y, rcond=None)[0]  # min ||y-P*x||_2
        r = y - np.dot(P, x)                      # update residual

    # embed coefficients in support
    x_hat = np.zeros((n, 1))
    x_hat[list(lamda)] = x

    return x_hat
