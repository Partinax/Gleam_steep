#!/usr/bin/env python

# Create a 'fake' TGSS image using the catalogue as input
# Convolve the fake image to the GLEAM resolution
# Subtract the fake image from the GLEAM postage stamp

from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
from scipy import signal
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from astropy import wcs
from optparse import OptionParser


usage="Usage: %prog [options] <file>\n"
parser = OptionParser(usage=usage)
parser.add_option('--catalogue',type="string", dest="catalogue",
                    help="The filename of the catalogue you want to read in.", default=None)
parser.add_option('--output',type="string", dest="output",
                    help="The name of the output catalogue.", default=None)
parser.add_option('--postage',dest="postage",default=None,type="string",
                  help="The filename of the postage stamp you want to read in.")
(options, args) = parser.parse_args()

cat = fits.open(options.catalogue)
cat_data = cat[1].data # treating an array as a method?

hdu = fits.open(options.postage)
header = hdu[0].header
w = wcs.WCS(header,naxis=2)# what does naxis=2 refer to?

try:
    delra = header["CD1_1"] # how is this exception handling working?
except:
    delra = header["CDELT1"] # why are these note referenced anywhere else?
try:
    deldec = header["CD2_2"]
except:
    deldec = header["CDELT2"]

naxis1 = header["NAXIS1"]
naxis2 = header["NAXIS2"]

print naxis1, naxis2 # what is naxis?

# I think wcs goes dec, ra:
#minra, mindec = w.all_pix2world(header["NAXIS2"],1,0)
#maxra, maxdec = w.all_pix2world(1,header["NAXIS1"],0)

# Slightly expand the selection because at the moment it's too small
#minra-=1
#maxra+=1
#mindec-=1
#maxdec+=1

#print minra, maxra, mindec, maxdec

#print minra, maxra, mindec, maxdec
coords = SkyCoord(cat_data["RA"], cat_data["Dec"], unit=(u.deg, u.deg)) 
x, y =  w.all_world2pix(coords.ra.deg,coords.dec.deg,0)
# boundaries are defined below
indices = np.where(x>0)
indices = np.intersect1d(indices,np.where(x<naxis1))
indices = np.intersect1d(indices,np.where(y>0))
indices = np.intersect1d(indices,np.where(y<naxis2))
indices = np.intersect1d(indices,np.where(np.isnan(cat_data["int_flux"])==False))
# Note to future self: select only sources which have a real int and peak flux (not nan or zero)

#if minra>180. and maxra<180.:
#    indices = np.where(cat_data["RA"]>maxra)
#    indices = np.concatenate((np.squeeze(indices),np.squeeze(np.where(cat_data["RA"]<minra))))
#else:
#    indices = np.where(cat_data["RA"]<maxra)
#    indices = np.intersect1d(indices,np.where(cat_data["RA"]>minra))
#indices = np.intersect1d(indices,np.where(cat_data["Dec"]<maxdec))
#indices = np.intersect1d(indices,np.where(cat_data["Dec"]>mindec))

cat[1].data = cat[1].data[indices]

cat.writeto(options.output,clobber=True)
