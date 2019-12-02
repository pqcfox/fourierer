from PIL import Image
import numpy as np
import argparse
import cmath

"""
Slow O(M^2N^2) implementation of the discrete Fourier transform.
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
    # for RGB images, which we don't have apparently
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
    m, n = imarray.shape
    result = np.zeros((m, n))
    
    if m % 2 != 0 or n % 2 != 0:
        raise ValueError('dimensions need to be a multiple of 2')
    elif m != n:
        raise ValueError('dimensions in this implementation need to be equivalent')
    elif m <= 32:
        return dft(imarray, inverse)
    else:
        # note that b/c of symmetry, the 2D DFT can be thought of as a 1D DFT horizontally, then a 1D DFT vertically
        for i in range(m):
            even = imarray[i][::2]
            odd = imarray[i][1::2]
            exponential = np.exp(-1j * 2 * np.pi * np.arange(n) / n)
            result[i] = np.concatenate([even + exponential[:n / 2] * odd, even + exponential[n / 2:] * odd])
        for j in range(n):
            even = result[:,j][::2]
            odd = imarray[:,j][1::2]
            exponential = np.exp(-1j * 2 * np.pi * np.arange(m) / m)
            result[:,j] = np.concatenate([even + exponential[:n / 2] * odd, even + exponential[n / 2:] * odd])
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
    m, n = image.size
    
    
def main():
    # load the data as a numpy array
    parser = argparse.ArgumentParser()
    parser.add_argument('image_file', type=str, help='path to the image file')
    parser.add_argument('implementation', type=str, help='which implementation of DFT to use, can be dft, c-t, or bruun')
    parser.add_argument('output_path', type=str, help='path to where the image should be saved')

    args = parser.parse_args()
    image = Image.open(args.image_file)
    imarray = np.array(image)

    if args.implementation == 'dft':
        # do a O(N^2) DFT on the image data
        convolved = dft(imarray)
    elif args.implementation == 'c-t':
        # use Cooley-Tukey divide-and-conquer approach
        convolved = ct(imarray)
    else:
        # use Bruun's implementation
        convolved = bruun(imarray)

    # save convolved image
    result = Image.fromarray(convolved)
    result.save('{}/output.png'.format(args.output_path), 'PNG')

if __name__ == '__main__':
    main()