# You are to implement the individual filter functions listed below. You can
# invoke this script in terminal using just one filter or multiple filters.
# Usage: python filter.py <path to input image> <filter name>\n"
# Multiple filters usage: python filter.py <path to input image> <filter 1> <filter 2> ..."
# Here are some examples:
# $ python filter.py keratocytes/A.tif mean
# 
# Or to use multiple filters, you can invoke it as such:
# $ python filter.py caulobacter/WT.tif threshold hp mean 
# The above applies a threshold filter to WT.tif, then applies a hp filter to
# the threshold output, and finally applies a mean filter to threshold-hp output. It
# is equivalent to mean(hp(threshold(WT.tif)))

from scipy import ndimage
from scipy.misc import imsave
import numpy as np
import matplotlib.pyplot as plt


def identity(X):
  return X

def mean(X):
  # TODO YOUR CODE HERE:
  # create an np.array R
  # return the result from ndimage.convolve(X,R)
  # where X is the original image
  # and R is the matrix which will be convolved with X

  R = np.array(()) # fill this in
  return None # replace with appropriate return value

  # END YOUR CODE HERE

def median(X):
  # TODO: observe the effect of changing the size of
  # the median filter
  size = 3
  return ndimage.median_filter(X,size)

def gaussian(X):
  # TODO: observe the effect of changing the std of
  # the gaussian filter
  std = 10
  return ndimage.gaussian_filter(X,std)

def hp(X):
  # TODO YOUR CODE HERE:
  # create an np.array R
  # return the result from ndimage.convolve(X,R)
  # where X is the original image
  # and R is the matrix which will be convolved with X
  R = np.array(())
  return None 
  # END YOUR CODE HERE

def gaussianHP(X):
  # TODO YOUR CODE HERE:
  # return an array that is the difference between
  # the original image and the Gaussian LP Filter
  return None
  # END YOUR CODE HERE

def threshold(X):
  # TODO YOUR CODE HERE:
  # Return a bit-mask where a pixel is 1 if it below
  # the threshold and 0 if it is above the threshold
  # (Note: you can write this in one line)
  return None
  # END YOUR CODE HERE


# This dictionary matches the name of each filter with its corresponding
# function (which you will implement above).
filters = {"identity":identity,
           "mean":mean,
           "median":median,
           "gaussian":gaussian,
           "threshold":threshold,
           "hp":hp,
           "gaussianHP":gaussianHP,
           }


def main():
  import sys,os
  usage = "\nUsage: python filter.py <path to input image> <filter name>\n"
  usage_multiple_filters = "\nMultiple filters usage: python filter.py <path to input image> <filter 1> <filter 2> ...\n"
  
  if len(sys.argv) < 3:
    print usage
    exit(1)

  # Checking valid path for image input file
  imgName = sys.argv[1]
  if not os.path.exists(imgName):
    print usage
    print "Error: File %s does not exist"%imgName
    exit(1)

  # Checking valid input for filter(s)
  filterNames = sys.argv[2:]
  for fn in filterNames:
    if fn not in filters:
      if len(sys.argv) == 3: print usage 
      else: print usage_multiple_filters
      print "Error: Filter '%s' was not found in filters\n" % f
      exit(1)

  # Reading in image
  img = plt.imread(imgName).astype(np.float64)
  if len(img.shape) > 2:
    print "Flattening image..."
    img = ndimage.imread(imgName, flatten=True)
  
  # Applying filter(s)
  fImg = img
  for fn in filterNames:
    f = filters[fn]
    fImg = f(fImg)
    # Checking if filter is implemented
    if fImg is None: 
      print "Your " + fn + " filter is not implemented."
      exit(1)
  
  # Printing image values to terminal
  print 'Original image: \n', img
  print 'Filtered image: \n', fImg

  # Using matplotlib to generate images colorbars
  plt.figure('original image')
  plt.imshow(img)
  plt.colorbar()
  plt.figure('filtered image')
  plt.imshow(fImg)
  plt.colorbar()
  plt.show()


# This boiler plate invokes the main() function when the script is run in
# python.
if __name__ == '__main__': 
  main()
