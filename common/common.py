#!/usr/bin/python

import matplotlib.pyplot as plt
import scipy.fftpack as fft
import numpy as np

def gen_test_signal(k=20, n=1000, amps_l=-100, amps_h=100, snr_db=None):
    (y, A, k, x_t, x_f) = gen_rand_gauss_signal(k, n, amps_l, amps_h, snr_db)
    return (y, A, x_t, x_f)

def add_awgn(x_volts, target_snr_db):
    # calculate signal power and convert to dB
    sig_avg_watts = np.mean(x_volts ** 2)
    sig_avg_db = 10 * np.log10(sig_avg_watts)

    # calculate noise according to SNR_dB = P_signal,dB - P_noise,dB then convert to watts
    noise_avg_db = sig_avg_db - target_snr_db
    noise_avg_watts = 10 ** (noise_avg_db / 10)

    # generate an sample of white noise
    mean_noise = 0
    noise_volts = np.random.normal(mean_noise, np.sqrt(noise_avg_watts), x_volts.shape)

    # noise up the original signal
    y_volts = x_volts + noise_volts

    return y_volts

def gen_rand_gauss_signal(k, n, amps_l, amps_h, snr_db):
    # number of samples needed
    c = 2
    m = int(np.ceil(c * k * np.log(n/k)))
    A = ((1.0/m)**0.5) * np.random.randn(m, n)

    # generate signal
    signal_f = np.zeros((n, 1))
    pos = np.random.randint(0, n, (k, 1))
    for i in range(0, len(pos)):
        signal_f[pos[i]] = np.random.randint(amps_l, amps_h)
    signal_t = fft.idct(signal_f, norm='ortho', axis=0)

    # add noise
    if snr_db:
        signal_t = add_awgn(signal_t, snr_db)
        signal_f = fft.dct(signal_t, norm='ortho', axis=0)

    plt_signal(signal_t, 'Time domain signal_t')
    plt_signal(signal_f, 'Frequency domain signal_f')

    # sample the signal
    y = np.dot(A, signal_t)

    D = fft.dct(np.eye(A.shape[1]), norm='ortho', axis=0)
    Ap = np.dot(A, D.T)

    return (y, Ap, k, signal_t, signal_f)

def plt_error(x_hat, x_f, title):
    # compute the error
    sse = np.sum((x_hat - x_f)**2)

    # plot it
    plt.figure(0)
    plt.stem(x_hat,  markerfmt='ro')
    plt.stem(x_f,  markerfmt='b-')
    plt.title(title + f', sse = {sse}')
    plt.show()

def plt_signal(signal, title):
    plt.plot(signal)
    plt.title(title)
    plt.grid()
    plt.show()
