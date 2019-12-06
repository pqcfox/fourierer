from PIL import Image
import numpy as np
import argparse
import cmath
import time

"""
Slow O(MN) implementation of the discrete Fourier transform.
Input:
    imarray = array representing the pixel intensities of 2D image
    inverse = boolean representing whether to take the inverse DFT
Output:
    result = array representing the resultant transformed image
"""
def dft(imarray, inverse=False):
    M, N = imarray.shape
    const = 1.0 / (M * N) if inverse else 1.0
    sign = 1.0 if inverse else -1.0

    result = np.zeros((M, N), dtype=complex)
    for k in range(M):
        for l in range(N):
            total = 0
            for a in range(M):
                for b in range(N):
                    multiplier = float(k * a) / M + float(l * b) / N
                    factor = cmath.exp(sign * 2.0j * cmath.pi * multiplier)
                    total += imarray[a][b] * factor * const
            result[k][l] = total
    return result


def row_column(image, fft_1d, inverse=False):
    M, N = image.shape
    output = np.zeros_like(image, dtype=complex)
    for i in range(N):
        output[i, :] = fft_1d(image[i, :], inverse=inverse)
    for j in range(N):
        output[:, j] = fft_1d(output[:, j], inverse=inverse)
    return output


def cooley_tukey_1d(X, inverse=False):
    N, = X.shape
    const = 1.0 / (M * N) if inverse else 1.0
    sign = 1.0 if inverse else -1.0

    if not np.log2(N).is_integer():
        raise ValueError("size of image must be a power of 2")

    if N == 1:
        return X

    X_even = cooley_tukey_1d(X[::2], inverse=inverse)
    X_odd = cooley_tukey_1d(X[1::2], inverse=inverse)
    twiddle = np.exp(-2.0j * np.pi * np.arange(N//2) / (N))
    X_odd_t = X_odd * twiddle
    plus = X_even + X_odd_t
    minus = X_even - X_odd_t
    result = np.concatenate((plus, minus))
    return result


def cooley_tukey(image, inverse=False):
    return row_column(image, cooley_tukey_1d, inverse=inverse)


def split_radix_1d(X, inverse=False):
    N, = X.shape
    const = 1.0 / (M * N) if inverse else 1.0
    sign = 1.0 if inverse else -1.0

    if not np.log2(N).is_integer():
        raise ValueError("size of image must be a power of 2")

    if N == 1:
        return X

    if N % 4 != 0:
        return cooley_tukey_1d(X, inverse=inverse)

    X_even = split_radix_1d(X[::2], inverse=inverse)
    X_one = split_radix_1d(X[1::4], inverse=inverse)
    X_three = split_radix_1d(X[3::4], inverse=inverse)

    X_even_front, X_even_back = X_even[:N//4], X_even[N//4:]
    twiddle_one = np.exp(2.0j * np.pi * np.arange(N//4) / N)
    twiddle_three = np.exp(-2.0j * np.pi * 3 * np.arange(N//4) / N)
    X_one_t = twiddle_one * X_one
    X_three_t = twiddle_three * X_three

    first = X_even_front + X_one_t + X_three_t
    second = X_even_back - 1.0j * (X_one_t - X_three_t)
    third = X_even_front - (X_one_t + X_three_t)
    fourth = X_even_back + 1.0j * (X_one_t - X_three_t)

    result = np.concatenate((first, second, third, fourth))
    return result


def split_radix(image, inverse=False):
    return row_column(image, split_radix_1d, inverse=inverse)


def vector_radix(X, inverse=False):
    M, N = X.shape
    const = 1.0 / (M * N) if inverse else 1.0
    sign = 1.0 if inverse else -1.0

    if M != N:
        raise ValueError("image must be square")

    if not np.log2(N).is_integer():
        raise ValueError("size of image must be a power of 2")

    if N == 1:
        return X

    X_even_even = vector_radix(X[::2, ::2], inverse=inverse)
    X_even_odd = vector_radix(X[::2, 1::2], inverse=inverse)
    X_odd_even = vector_radix(X[1::2, ::2], inverse=inverse)
    X_odd_odd = vector_radix(X[1::2, 1::2], inverse=inverse)
    f_even_odd = np.exp(sign * 2.0j * np.pi * np.arange(N) / N)[None, :]
    f_odd_even = np.exp(sign * 2.0j * np.pi * np.arange(N) / N)[:, None]
    f_odd_odd = f_even_odd * f_odd_even

    result = np.tile(X_even_even, (2, 2)).astype(complex)
    result += f_even_odd * np.tile(X_even_odd, (2, 2))
    result += f_odd_even * np.tile(X_odd_even, (2, 2))
    result += f_odd_odd * np.tile(X_odd_odd, (2, 2))
    return result


def pfa_1d(x):
    pass


"""
Popular divide-and-conquer approach to DFT proposed by Cooley and Tukey.
Implements a radix-2 decimation-in-time row-column Cooley-Tukey FFT.
Input:
    imarray = array representing the pixel intensities of 2D image
Output:
    result = array representing the DFT of the image
"""
def ct_fast(imarray, inverse=False):
    N = imarray.shape[0]
    const = 1.0 / (M * N) if inverse else 1.0
    sign = 1.0 if inverse else -1.0

    if not np.log2(N).is_integer():
        raise ValueError("size of image must be a power of 2")

    # N_min here is equivalent to the stopping condition above,
    # and should be a power of 2
    N_min = min(N, 32)

    def ct_1d(x):
        # Perform an O[N^2] DFT on all length-N_min sub-problems at once
        k = np.arange(N_min)[None, :]
        n = np.arange(N_min)[:, None]
        M = np.exp(sign * 2.0j * np.pi * n * k / N_min)
        X = np.dot(M, x.reshape((N_min, -1)))

        # build-up each level of the recursive calculation all at once
        while X.shape[0] < N:
            X_even = X[:, :(int(X.shape[1] / 2))]
            X_odd = X[:, (int(X.shape[1] / 2)):]
            factor = np.exp(sign * 1.0j * np.pi * np.arange(X.shape[0])
                            / X.shape[0])[:, None]
            X = np.vstack([X_even + factor * X_odd,
                        X_even - factor * X_odd])

        return X.ravel()

    inter = np.zeros((N, N), dtype=complex)
    for n in range(N):
        x = imarray[n]
        inter[n] = ct_1d(x)
    result = np.zeros((N, N), dtype=complex)
    for n in range(N):
        x = inter[:,n]
        result[:,n] = ct_1d(x)

    return result




def main():
    # load the data as a numpy array
    parser = argparse.ArgumentParser()
    parser.add_argument('image_file', type=str, help='path to the image file')
    parser.add_argument('implementation', type=str, help='which implementation of DFT to use, can be dft, c-t, or bruun')
    parser.add_argument('output_path', type=str, help='path to where the image should be saved')

    args = parser.parse_args()
    image = Image.open(args.image_file)
    imarray = np.array(image)

    result = None
    start_time = time.time()
    if args.implementation == 'dft':
        # do a O(N^4) DFT on the image data
        result = dft(imarray)
    elif args.implementation == 'ct':
        # use Cooley-Tukey divide-and-conquer approach
        result = cooley_tukey(imarray)
    elif args.implementation == 'sr':
        # use split-radix implementation
        result = split_radix(imarray)
    elif args.implementation == 'vr':
        # use vector-radix implementation
        result = vector_radix(imarray)
    else:
        raise ValueError('not a valid implementation')
    print("--- %s seconds ---" % (time.time() - start_time))

    # to display the image, we only take the magnitude of the resultant array, as it contains most of the geometric structure info
    magnitude = np.real(result)

    # save convolved image
    im_result = Image.fromarray(magnitude).convert('L')
    im_result.save(args.output_path, 'PNG')

    # get numpy fft result
    start_time = time.time()
    result = np.fft.fft2(imarray)
    print("--- %s seconds ---" % (time.time() - start_time))

    magnitude = np.real(result)
    im_result = Image.fromarray(magnitude).convert('L')
    im_result.save('true.png', 'PNG')



if __name__ == '__main__':
    main()
