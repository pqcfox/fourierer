#  This is a simple script to generate images using matplotlib, a python
#  plotting library. To view, for example, one of the keratocyte images you
#  would invoke the script as follows in the terminal: 
#  python show.py keratocytes/A.tif  

#  This script also prints the raw values of the pixel values comprising the
#  image (in the format of a 2D array)

import matplotlib.pyplot as plt
import sys

usage = "Usage: python show.py <path to image file>"

if len(sys.argv) < 1:
  print usage
  exit(1)

img = plt.imread(sys.argv[1])
print 'Image values: '
print img # This will print the raw pixel values for the image

plt.rcParams['image.cmap']='jet'
plt.imshow(img)
plt.colorbar()
plt.show()
