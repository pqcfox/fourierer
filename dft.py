from PIL import Image
import numpy as np
import argparse
import cmath

"""
Slow O(M^2N^2) implementation of the discrete Fourier transform.
Input:
    image = PIL Image object
    inverse = boolean representing whether to take the inverse DFT
Output:
    dft_r, dft_g, dft_b = matrices representing pixel intensities at each red, green, and blue channel
"""
def dft(image, inverse=False):
    const = 1
    if inverse == True:
        const = 1 / (2 * pi)
    M, N = image.size
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
    result = Image.fromarray(numpy.T, 'RGB')
    return result


"""
Popular divide-and-conquer approach to DFT proposed by Cooley and Tukey.
Input:
    image = PIL Image object
Output:
    dft_r, dft_g, dft_b = matrices representing pixel intensities at each red, green, and blue channel
"""
def ct(image):
     M, N = image.size


"""
Bruun's FFT using a recursive polynomial-factorization approach. Must be used on factors of 2.
Input:
    image = PIL Image object
Output:
    dft_r, dft_g, dft_b = matrices representing pixel intensities at each red, green, and blue channel
"""
def bruun(image):
    
def main():
    # load the data as a numpy array
    parser = argparse.ArgumentParser()
    parser.add_argument('image_file', type=str, help='path to the image file')
    parser.add_argument('implementation', type=str, help='which implementation of DFT to use, can be dft, c-t, or bruun')
    parser.add_argument('output_path', type=str, help='path to where the image should be saved')

    args = parser.parse_args()
    image = Image.open(args.image_file)

    if args.implementation == 'dft':
        # do a O(N^2) DFT on the image data
        convolved = dft(image)
    elif args.implementation == 'c-t':
        # use Cooley-Tukey divide-and-conquer approach
        convolved = ct(image)
    else:
        # use Bruun's implementation
        convolved = bruun(image)

    # save convolved image
    convolved.save('{}/output.png'.format(args.output_path), 'PNG')

if __name__ == '__main__':
    main()