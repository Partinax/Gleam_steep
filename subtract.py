#!/usr/bin/env python

# Subtract one fits file from another

import sys

# fits files
try:
    import astropy.io.fits as fits
except ImportError:
    import pyfits as fits

file1=sys.argv[1]
file2=sys.argv[2]
output=sys.argv[3]

print "Subtracting "+file2+" from "+file1

hdu1=fits.open(file1)
hdu2=fits.open(file2)
hdu1[0].data-=hdu2[0].data
hdu1.writeto(output,clobber=True)
