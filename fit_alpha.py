#!/usr/bin/env python

# Fit simple power-law SEDs to the sources detected in the red image

import os, sys

# Need a least-squares estimator that gives a useable error estimate
from scipy.optimize import leastsq

import numpy as np
#tables and votables
import astropy.io.fits as fits
from astropy.io.votable import parse_single_table
from astropy.io.votable import writeto as writetoVO
from astropy.table import Table, Column

# Parallelise the code
import multiprocessing
# multiple cores support
import pprocess

import matplotlib.pyplot as plt

from optparse import OptionParser

usage="Usage: %prog [options] <file>\n"
parser = OptionParser(usage=usage)
parser.add_option('--catalogue',type="string", dest="catalogue",
                    help="The filename of the catalogue you want to read in.", default=None)
parser.add_option('--output',type="string", dest="output",
                    help="The filename of the output.", default=None)
#parser.add_option('--andre',type="string", dest="andre",
#                    help="If desired, the output text file suitable for use in Andre's tools (default = don't make one).", default=None)
#parser.add_option('--plot',action="store_true",dest="make_plots",default=False,
#                  help="Make fit plots? (default = False)")
#parser.add_option('--order',dest="poly_order",default=1,type=int,
#                  help="Set the order of the polynomial fit. (default = 1)")
parser.add_option('--cores',dest="cores",default=None,type=int,
                  help="How many cores to use. (default = all available)")
parser.add_option('--limit',dest="flux_limit",default=0.0,type=float,
                  help="Minimum flux density at any frequency for source to be included to in fit (default = 0.0Jy)")
parser.add_option('--calerror',dest="calerror",default=0.02,type=float,
                  help="Estimated calibration error for -72<Dec<+18.5 (default = 0.02, i.e. 2%)")
parser.add_option('--hldec_calerror',dest="hldec_calerror",default=0.03,type=float,
                  help="Estimated calibration error for Dec<-72 and Dec>18.5 (default = 0.03, i.e. 3%)")
(options, args) = parser.parse_args()

# http://scipy-cookbook.readthedocs.org/items/FittingData.html
# Define function for calculating a power law
def powerlaw(x,amp,index):
    return amp * (x**index)
vpowerlaw = np.vectorize(powerlaw)

# define our (line) fitting function
fitfunc = lambda p, x: p[0] + p[1] * x
errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

if options.cores is None:
    cores = multiprocessing.cpu_count()
else:
    cores = options.cores
if options.output is None:
    output="test.vot"
else:
    output=options.output

if options.catalogue is None:
    print "must specify input catalogue"
    sys.exit(1)
else:
    filename, file_extension = os.path.splitext(options.catalogue)
    if file_extension == ".fits":
        temp = fits.open(options.catalogue)
        data = temp[1].data
        temp.close()
    elif file_extension == ".vot":
        temp = parse_single_table(options.catalogue)
        data = temp.array

# Spectral index fitting function

def fit_spectrum(freq_array,flux_array,flux_errors): #,plot):
    pinit = [-2.0, -0.7]
    pinit = [0.0, -0.7]
    fit = leastsq(errfunc, pinit, args=(freq_array, flux_array, flux_errors), full_output=1)
    covar = fit[1]
    if covar is not None:
        P = fit[0]
        residual = errfunc(P,freq_array, flux_array, flux_errors)
        chi2red = sum(np.power(residual,2))/(len(freqs)-len(pinit))
        alpha=P[1]
        amp = np.exp(P[0])
    # Errors
        err_alpha = np.sqrt(covar[0][0])
        err_amp = np.sqrt(covar[1][1])
    else:
        chi2red=None
        alpha=None
        amp=None
        err_alpha=None
        err_amp=None
#    if plot:
#        outpng=(name.replace(" ","_"))+".png"
#        # Plot
#        example=plt.figure(figsize=(10,5))
#        example.suptitle(name)
#        ax1=example.add_subplot(1,2,1)
#        ax1.plot(np.exp(freq_array), powerlaw(np.exp(freq_array), amp, alpha),label="alpha={0:3.1f} ".format(alpha))     # Scipy Fit
#        ax1.errorbar(np.exp(freq_array), np.exp(flux_array), yerr=flux_errors*np.exp(flux_array), fmt='k.')  # Data
#        ax1.set_ylabel("S / Jy")
#        ax1.set_xlabel("Frequency / MHz")
#        ax1.legend()
#        ax2=example.add_subplot(1,2,2)
#        ax2.loglog(np.exp(freq_array), powerlaw(np.exp(freq_array), amp, alpha),label="alpha={0:3.1f} ".format(alpha))     # Scipy Fit
#        ax2.errorbar(np.exp(freq_array), np.exp(flux_array), yerr=flux_errors*np.exp(flux_array), fmt='k.')  # Data
##        ax2.set_xlim(left=xmin01,right=xmax01)
##       ax2.set_ylim([min(np.exp(flux_array)),max(np.exp(flux_array))])
#        ax2.set_ylabel("S / Jy")
#        ax2.set_xlabel("Frequency / MHz")
#        ax2.legend()
#        example.savefig(outpng)
#        example.clf()
    return alpha, err_alpha, amp, err_amp, chi2red

# Set up the parameters and do the fitting

# Frequencies to write out -- use the full band since it shows the range I did the fitting over
freq1=72
freq2=231

# Can't figure out how to read the bloody column names! hardcode
#freqs=["076", "084", "092", "099", "107",  "115", "122", "130", "143", "151", "158", "166", "174",  "181", "189","197", "204", "212", "220", "227"]
freqs=["088","118","154","200"]
# Frequencies
freq_array = np.log(np.ma.array([float(int(x)) for x in freqs]),dtype="float32")

# Select source flux densities that meet the flux limit criterion
#if options.flux_limit is not None:
#    brightsrcs = np.where(data["int_flux_wide"]>options.flux_limit)
#else:
#    brightsrcs = np.where(data["int_flux_wide"]>0.0)
# And are not NaN at any frequency
brightsrcs = np.where(data["int_flux_088"]>options.flux_limit)
#for x in freqs:
#    brightsrcs = np.intersect1d(brightsrcs,np.where(np.isfinite(data["int_flux_"+x])))
#    brightsrcs = np.intersect1d(brightsrcs,np.where(data["int_flux_"+x]>options.flux_limit))

brightsrcs = np.squeeze(brightsrcs)

print "Fitting",len(brightsrcs),"source spectral energy distributions"

# Extreme Dec sources = 3% calibration error
calibration_error = options.hldec_calerror*np.ones(len(brightsrcs))

# Mid-range Dec sources = 2% calibration error
midsrcs = np.intersect1d(brightsrcs,np.where(data["dec_088"]<=18.5))
midsrcs = np.intersect1d(midsrcs,np.where(data["dec_088"]>=-72.))
midsrcs = np.squeeze(midsrcs)
calibration_error[midsrcs] = options.calerror

flux_list = []
err_list = []
for x in freqs:
    flux_list.append(data["int_flux_"+x][brightsrcs])
# Error propagation: error on log(x) = err_x/x
    fitting_error = data["err_int_flux_"+x][brightsrcs]/data["int_flux_"+x][brightsrcs]
    err_list.append(np.sqrt(fitting_error**2 + calibration_error**2))
    
flux_array = np.transpose(np.ma.vstack(flux_list)).astype("float32")
flux_array = np.ma.log(flux_array)
flux_errors = np.transpose(np.ma.vstack(err_list)).astype("float32")
#names = data["Name"][brightsrcs]

#weights = 1/(flux_errors*flux_errors)

results = pprocess.Map(limit=cores)
calc = results.manage(pprocess.MakeParallel(fit_spectrum))

for i in range(0,len(brightsrcs)):
    calc(freq_array,flux_array[i],flux_errors[i]) # ,options.plot)

# Unpack results
alpha, err_alpha, amp, err_amp, chi2red = map(list, zip(*results))

# Convert to numpy arrays
alpha = np.array(alpha, dtype="float32")
err_alpha = np.array(err_alpha, dtype="float32")
amp = np.array(amp, dtype="float32")
err_amp = np.array(err_amp, dtype="float32")
chi2red = np.array(chi2red, dtype="float32")

# Exclude any sources which came out with NaN alphas or amps
#good = np.squeeze(np.where(np.isfinite(alpha)))

# Generate flux density columns
#flux1 = vpowerlaw(np.tile(freq1,len(amp[good])),amp[good],alpha[good])
#err_flux1 = err_amp[good]*flux1
#flux2 = vpowerlaw(np.tile(freq2,len(amp[good])),amp[good],alpha[good])
#err_flux2 = err_amp[good]*flux2

# Add a column to the FITS table
table = Table.read(options.catalogue)
alpha_col = Column(name='alpha', data=alpha)
err_alpha_col = Column(name='err_alpha', data=err_alpha)
S88_col = Column(name='int_flux_fit_88', data=vpowerlaw(np.tile(88.,len(amp)),amp,alpha))
err_S88_col = Column(name='err_int_flux_fit_88', data=err_amp*88.)
chi2red_col = Column(name='chi2red', data=chi2red)
table.add_column(alpha_col)
table.add_column(err_alpha_col)
table.add_column(S88_col)
table.add_column(err_S88_col)
table.add_column(chi2red_col)
table.write("new.fits",format="fits")
