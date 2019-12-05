from PIL import Image
import numpy as np
import argparse
import cmath

"""
Slow O(MN) implementation of the discrete Fourier transform.
Input:
    imarray = array representing the pixel intensities of 2D image
    inverse = boolean representing whether to take the inverse DFT
Output:
    dft_r, dft_g, dft_b = matrices representing pixel intensities at each red, green, and blue channel
"""
def dft(imarray, inverse=False):
    const = 1
    if inverse == True:
        const = 1 / (2 * pi)
    m, n = imarray.shape
    return np.array([[sum([sum([imarray[i,j] * np.exp(-1j*2*np.pi*(k_m*i/m + k_n*j/n)) * const for i in range(m)]) for j in range(n)]) for k_n in range(n)] for k_m in range(m)])

    """
    # for RGB images
    dft_r = np.zeros((M, N))
    dft_g = np.zeros((M, N))
    dft_b = np.zeros((M, N))
    pixels = image.load()
    for i in range(M):
        for j in range(N):
            sum_r = 0
            sum_g = 0
            sum_b = 0
            for m in range(M):
                for n in range(N):
                    (r, g, b, alpha) = pixels[m,n]
                    exponential = cmath.exp(-2 * cmath.pi * 1j * (float(k * m) / M + float(l * n) / N))
                    sum_r += r * exponential * const
                    sum_g += g * exponential * const
                    sum_b += b * exponential * const
            dft_r[l][k] = sum_r / M / N
            dft_g[l][k] = sum_g / M / N
            dft_b[l][k] = sum_b / M / N
    result = np.array([dft_r, dft_g, dft_b])
    print(result.shape)
    result = Image.fromarray(np.T, 'RGB')
    return result
    """


"""
Popular divide-and-conquer approach to DFT proposed by Cooley and Tukey.
Input:
    imarray = array representing the pixel intensities of 2D image
    inverse = boolean representing whether to take the inverse DFT
Output:
    result = array representing the convolved image
"""
def ct(imarray, inverse=False):
    const = 1
    if inverse == True:
        const = 1 / (2 * pi)
    N = imarray.shape[0]
    
    if np.log2(N) % 1 > 0:
        raise ValueError("size of image must be a power of 2")

    # N_min here is equivalent to the stopping condition above,
    # and should be a power of 2
    N_min = min(N, 32)

    def ct_1d(x):
        # Perform an O[N^2] DFT on all length-N_min sub-problems at once
        n = np.arange(N_min)
        k = n[:, None]
        M = np.exp(-2j * np.pi * n * k / N_min)
        X = np.dot(M, x.reshape((N_min, -1)))

        # build-up each level of the recursive calculation all at once
        while X.shape[0] < N:
            X_even = X[:, :(int(X.shape[1] / 2))]
            X_odd = X[:, (int(X.shape[1] / 2)):]
            factor = np.exp(-1j * np.pi * np.arange(X.shape[0])
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


"""
Bruun's FFT using a recursive polynomial-factorization approach. Must be used on factors of 2.
Input:
    imarray = array representing the pixel intensities of 2D image
    inverse = boolean representing whether to take the inverse DFT
Output:
    result = array representing the convolved image
"""
def bruun(imarray, inverse=False):
    m, n = imarray.size
    return None
    
    
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
    if args.implementation == 'dft':
        # do a O(N^4) DFT on the image data
        result = dft(imarray)
    elif args.implementation == 'c-t':
        # use Cooley-Tukey divide-and-conquer approach
        result = ct(imarray)
    elif args.implementation == 'bruun':
        # use Bruun's implementation
        result = bruun(imarray)
    else:
        raise ValueError('not a valid implementation')

    # to display the image, we only take the magnitude of the resultant array, as it contains most of the geometric structure info
    magnitude = np.real(result)

    # save convolved image
    im_result = Image.fromarray(magnitude).convert('L')
    im_result.save(args.output_path, 'PNG')

    # get numpy fft result
    result = np.fft.fft2(imarray)
    magnitude = np.real(result)
    im_result = Image.fromarray(magnitude).convert('L')
    im_result.save('true.png', 'PNG') 



if __name__ == '__main__':
    main()