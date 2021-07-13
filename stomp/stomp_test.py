#!/usr/bin/python

import matplotlib.pyplot as plt
import scipy.fftpack as fft
import numpy as np

from stomp import *

def gen_rand_gauss_signal(k=20, n=10000, amps_l=-100, amps_h=100):
    # number of samples needed
    c = 4
    m = int(np.ceil(c * k * np.log(n/k)))
    A = ((1.0/m)**0.5) * np.random.randn(m, n)

    # generate signal
    signal_f = np.zeros((n, 1))
    pos = np.random.randint(1, n, (k, 1))
    for i in range(0, len(pos)):
        signal_f[pos[i]] = np.random.randint(amps_l, amps_h)
    signal_t = fft.idct(signal_f, norm='ortho', axis=0)

    # sample the signal
    y = np.dot(A, signal_t)

    return (y, A, k, signal_t, signal_f)

def plt_error(x_f, x_hat):
    # compute the error
    sse = np.sum((x_hat - x_f)**2)
    print(f"sse = {sse}")

    # plot it
    plt.figure(0)
    plt.stem(x_hat,  markerfmt='o')
    plt.stem(x_f,  markerfmt='-')

    plt.figure(1)
    plt.stem(x_hat - x_f)
    plt.show()

def stomp_test():
    sigma = 3  # check for 3 std devs away
    stages = 3  # n stages in stomp

    # -------------------------------------------------------------------------------- 
    (y, A, k, x_t, x_f) = gen_rand_gauss_signal() # make signal

    print(f"calling stomp()")

    D = fft.dct(np.eye(A.shape[1]), norm='ortho', axis=0)      # optimize over the dct domain
    x_hat = stomp(y, np.dot(A, D.T), sigma, stages)            # stomp call with n stages
    plt_error(x_f, x_hat)                                      # print and graph errors

if __name__ == "__main__":
    stomp_test()