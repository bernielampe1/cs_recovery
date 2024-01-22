#!/usr/bin/python

"""

References:
"""

import numpy as np
import scipy.linalg as lin

def adamm(y, A, term, param, max_iters):
    """
    Parameters:
        y: `compressed samples`
        A: `sampling matrix`
        term: `termination function`
        param: `termination parameter` { k: `sparsity` | p: `percent` | e: `epsilon` }

    Returns:
        x_hat: `reconstructed signal`
    """

    k = 0
    N = A.shape[1]
    x_k = np.zeros((N, 1))
    v_k = np.zeros((N, 1))
    u_k = np.zeros((N, 1))
    I = np.identify(N)

    while not term(param, y, r, x_hat) and k < max_iters:
        alpha = np.inv(np.dot(A.T, A + I))
        beta = np.dot(A.T, y) + v_k - u_k
        x_k = np.dot(alpha, beta)

        v_k = S_delta(x_k + u_k)

        u_k = u_k + x_k - v_k

    return x_k

# terminate with output signal has sparsity k
def sparsity_term(k, y, r, x_hat):
    return int(np.linalg.norm(x_hat, 0, axis=0)) == k

# terminate when output signal has p percentage of signal
def percent_term(p, y, r, x_hat):
    y_l2p = np.linalg.norm(y) * (1.0-p)
    return y_l2p > np.linalg.norm(r)

# terminate when energy (L2 norm) of remaining residual falls below this energy
def epsilon_term(e, y, r, x_hat):
    return e > np.linalg.norm(r)

if __name__ == "__main__":
    from common import *

    max_iters = 100

    for db in [None, 20]:
        (y, A, x_t, x_f) = gen_test_signal(snr_db=db)

        x_h = mp(y, A, sparsity_term, 20, max_iters)
        plt_error(x_h, x_f, 'sparsity_term, k = 20', max_iters)

        x_h = mp(y, A, percent_term, 0.99)
        plt_error(x_h, x_f, 'percent_term, p = 0.99', max_iters)

        x_h = mp(y, A, epsilon_term, 0.00001)
        plt_error(x_h, x_f, 'epsilon_term, k = 0.00001', max_iters)
