#!/usr/bin/python

# ref0: https://en.wikipedia.org/wiki/Conjugate_gradient_method
# ref1: l1magic, Justin Romberg
# reg2: Shewchuk, Jonathan R (1994). An Introduction to the Conjugate Gradient Method Without the Agonizing Pain

import numpy as np

def cgsolve(A, b, tol = 1e-8, maxiter = 1000, verbose = 0):
    x = np.zeros(b.shape)
    r = b # - A @ x, but x_0 = vec(0)
    d = r
    delta = np.dot(r.T, r)
    delta0 = np.dot(b.T, b)
    numiter = 0
    bestx = x
    bestres = np.sqrt(delta/delta0)
    while numiter < maxiter and delta > tol**2 * delta0:

        q = A @ d
        alpha = delta / np.dot(d.T, q)
        x = x + alpha * d

        if np.mod(numiter, 50) == 0: # from jonathan shewchuk, page 50, section B2. Conjugate Gradients
            r = b - A @ x
        else:
            r = r - alpha * q

        deltaold = delta
        delta = np.dot(r.T, r)
        beta = delta / deltaold
        d = r + beta * d
        numiter += 1

        curres = np.sqrt(delta/delta0) # from l1-magic by justin romberg
        if curres < bestres: # keep best answer so far incase of error due to round off accumulation
            bestx = x
            bestres = curres

        if verbose:
            print(f'cg: iter = {numiter}, best res = {bestres}, cur res = {curres}')
    return bestx

if __name__ == "__main__":
    n = 100 # test dimension
    s = 10 # test scale
    A = np.random.rand(n, n) * s # scaled random matrix
    A = A @ A.T # make psd
    b = np.random.rand(n) * s
    x = cgsolve(A, b, verbose=1)
    r = np.linalg.norm(A @ x - b)
    print(f'residual L2 norm = {r}')

