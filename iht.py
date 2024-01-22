#!/usr/bin/python

""" Iterative Hard Thresholding for Compressive Sensing
https://arxiv.org/pdf/0805.0510.pdf

References:
    Blumensath, Thomas, and Mike E. Davies. "Iterative hard thresholding for compressed sensing."
    Applied and computational harmonic analysis 27.3 (2009): 265-274.
"""

import numpy as np
import scipy.linalg as lin

def iht(y, A, k, tol):
    """
    Parameters:
        y: `compressed samples`
        A: `sampling matrix`
        k: `sparity`
        tol: `smallest residual to stop`

    Returns:
        x_hat: `reconstructed signal`
    """

    r = np.copy(y)              # init residual
    N = A.shape[1]
    x_hat = np.zeros((N, 1))    # null output init

    i = 0
    while i < 10:
        # compute correlations
        beta = x_hat + np.dot(A.T, r)

        # hard threshold
        inds = np.abs(beta) < np.partition(beta.flatten(), -k)[-k]
        beta[inds] = 0
        x_hat = np.copy(beta)

        # update residual
        r = y - np.dot(A, x_hat)

        i += 1

    return x_hat

if __name__ == "__main__":
    from common import *

    for db in [None, 20]:
        (y, A, x_t, x_f) = gen_test_signal(snr_db=db)

        x_h = iht(y, A, 20, 0.0001)
        plt_error(x_h, x_f, 'sparsity_term, k = 20')
