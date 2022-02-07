"""
Compressive Sampling Matching Pursuit

References:


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
        sigma: `std deviation threshold (2-3)`
        term: `termination criteria function`
        param: `termination criteria parameter`

    Returns:
        x_hat: `reconstructed signal`
    """

    r = y       # init residual
    lamda = set()  # list of support loc inds
    x = np.array([]) # null output init
    n = A.shape[1] # size of solution

    for t in range(N):
        c = np.dot(A.T, r)                  # correlation
        thr = sigma * np.std(r)             # compute number of std devs for threshold
        inds = np.nonzero(np.abs(c) > thr)  # find abs support above threshold
        lamda.update(inds[0].tolist())      # update support

        P = A[:, list(lamda)]               # compose the submatrix
        x = np.linalg.lstsq(P, y, rcond=None)[0]  # min ||y-P*x||_2
        r = y - np.dot(P, x)                      # update residual

    # embed coefficients in support
    x_hat = np.zeros((n, 1))
    x_hat[list(lamda)] = x

    return x_hat
