#!/usr/bin/python

""" Basic FOCUSS algorithm

References:
    I. F. Gorodnitsky and B. D. Rao, "Sparse signal reconstruction from limited data using FOCUSS: a re-weighted minimum norm algorithm," in IEEE Transactions on Signal Processing, vol. 45, no. 3, pp. 600-616, March 1997, doi: 10.1109/78.558475.
"""

import numpy as np

def focuss(y, A, itrs):
    """
    Parameters:
        y: `compressed samples`
        A: `sampling matrix`
        itrs: `number of algorithm itrs`

    Returns:
        x_hat: `reconstructed signal`
    """

    x_k = np.ones((A.shape[1], 1))  # init solution non-zero

    i = 0
    while i < itrs:
        W_pk = np.diagflat(x_k)             # Step 1: W_pk = diag(x_k-1)
        alpha = np.dot(A, W_pk)             # Step 2: q_k = (A*W_pk)^+y
        alpha_plus = np.linalg.pinv(alpha) 
        q_k = np.dot(alpha_plus, y)
        x_k = np.dot(W_pk, q_k)             # Step 3: x_k = W_pk * q_k

        i+=1

    return x_k


if __name__ == "__main__":
    from common import *

    for db in [None, 20]:
        (y, A, x_t, x_f) = gen_test_signal(snr_db=db)

        x_h = focuss(y, A, 500)
        plt_error(x_h, x_f, 'sparsity_term, k = 20')

