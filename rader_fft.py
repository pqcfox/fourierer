import numpy as np


def rader_1d(X):
    N, = X.shape

    # Ensure that our input is prime
    if any([N % i == 0 for i in range(2, int(np.sqrt(N)) + 1)]):
        raise ValueError("size of input must be prime")

    # Find a primitive root mod N
    for g in range(2, N):
        for order in range(1, N):
            if pow(g, order, N) == 1:
                break

        if order == N - 1:
            break

    # Set the first entry to the sum of the inputs
    result = np.zeros_like(X, dtype=complex)
    result[0] = np.sum(X)

    # Calculate the exponents of g mod q
    exps = np.array([pow(g, q, N) for q in range(N - 1)])
    neg_exps = np.array([pow(g, q, N) for q in range(N - 1, 0, -1)])

    # Find the two sequences we convolve
    a = X[neg_exps]
    b = np.exp(-2.0j * np.pi * exps / N)

    # Convolve them with the built-in numpy FFT
    # (this can be any FFT algorithm)
    prod = np.fft.fft(b) * np.fft.fft(a)
    conv = np.fft.ifft(prod)

    # Combine results and return
    result[exps] = X[0] + conv
    return result


X = np.array([3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9])
print('Rader\'s FFT: {}'.format(rader_1d(X)))
print('Numpy FFT: {}'.format(np.fft.fft(X)))
