#!/bin/bash

input=~/Dropbox/Public/GLEAM_EGC.fits
output=./GLEAM_EGC_aegean_columns.fits

stilts tpipe in=$input cmd='addcol ra "RAJ2000"' \
                       cmd='addcol err_ra "err_RAJ2000"' \
                       cmd='addcol dec "DEJ2000"' \
                       cmd='addcol err_dec "err_DEJ2000"' \
                       cmd='addcol int_flux "(int_flux_076 + int_flux_084 + int_flux_092 + int_flux_099)/4.0"' \
                       cmd='addcol err_int_flux "sqrt(err_int_flux_076*err_int_flux_076 + err_int_flux_084*err_int_flux_084 + err_int_flux_092*err_int_flux_092 + err_int_flux_099*err_int_flux_099)"' \
                       cmd='addcol peak_flux "(peak_flux_076 + peak_flux_084 + peak_flux_092 + peak_flux_099)/4.0"' \
                       cmd='addcol err_peak_flux "sqrt(err_peak_flux_076*err_peak_flux_076 + err_peak_flux_084*err_peak_flux_084 + err_peak_flux_092*err_peak_flux_092 + err_peak_flux_099*err_peak_flux_099)"' \
                       cmd='addcol a "a_076"' \
                       cmd='addcol err_a "a_076*err_a_wide/a_wide"' \
                       cmd='addcol b "b_076"' \
                       cmd='addcol err_b "b_076*err_b_wide/b_wide"' \
                       cmd='addcol pa "pa_076"' \
                       cmd='addcol err_pa "pa_076*err_pa_wide/pa_wide"' \
                       cmd='addcol local_rms "sqrt(local_rms_076*local_rms_076 + local_rms_084*local_rms_084 + local_rms_092*local_rms_092 + local_rms_099*local_rms_099)"' \
                       cmd='keepcols "ra err_ra dec err_dec int_flux err_int_flux peak_flux err_peak_flux a err_a b err_b pa err_pa local_rms"' \
                       out=$output

#                       out=temp1.fits
#
#stilts tpipe in=$input cmd='select !NULL_alpha' \
#                       cmd='addcol ra "RAJ2000"' \
#                       cmd='addcol err_ra "err_RAJ2000"' \
#                       cmd='addcol dec "DEJ2000"' \
#                       cmd='addcol err_dec "err_DEJ2000"' \
#                       cmd='addcol int_flux "(int_flux_076*pow(88./76.,-0.83) + int_flux_084*pow(88./84.,-0.83) + int_flux_092*pow(88./92.,-0.83) + int_flux_099*pow(88./99.,-0.83))/4.0"' \
#                       cmd='addcol err_int_flux "sqrt(err_int_flux_076*err_int_flux_076 + err_int_flux_084*err_int_flux_084 + err_int_flux_092*err_int_flux_092 + err_int_flux_099*err_int_flux_099)"' \
#                       cmd='addcol peak_flux "(peak_flux_076*pow(88./76.,-0.83) + peak_flux_084*pow(88./84.,-0.83) + peak_flux_092*pow(88./92.,-0.83) + peak_flux_099*pow(88./99.,-0.83))/4.0"' \
#                       cmd='addcol err_peak_flux "sqrt(err_peak_flux_076*err_peak_flux_076 + err_peak_flux_084*err_peak_flux_084 + err_peak_flux_092*err_peak_flux_092 + err_peak_flux_099*err_peak_flux_099)"' \
#                       cmd='addcol a "a_076"' \
#                       cmd='addcol err_a "a_076*err_a_wide/a_wide"' \
#                       cmd='addcol b "b_076"' \
#                       cmd='addcol err_b "b_076*err_b_wide/b_wide"' \
#                       cmd='addcol pa "pa_076"' \
#                       cmd='addcol err_pa "pa_076*err_pa_wide/pa_wide"' \
#                       cmd='keepcols "ra err_ra dec err_dec int_flux err_int_flux peak_flux err_peak_flux a err_a b err_b pa err_pa"' \
#                       out=temp2.fits

#stilts tcat in='temp1.fits temp2.fits' out=$output

#                       cmd='addcol int_flux "(int_flux_076*pow(88./76.,alpha) + int_flux_084*pow(88./84.,alpha) + int_flux_092*pow(88./92.,alpha) + int_flux_099*pow(88./99.,alpha))/4.0"' \
#                       cmd='addcol err_int_flux "sqrt(err_int_flux_076*err_int_flux_076 + err_int_flux_084*err_int_flux_084 + err_int_flux_092*err_int_flux_092 + err_int_flux_099*err_int_flux_099)"' \
#                       cmd='addcol peak_flux "(peak_flux_076*pow(88./76.,alpha) + peak_flux_084*pow(88./84.,alpha) + peak_flux_092*pow(88./92.,alpha) + peak_flux_099*pow(88./99.,alpha))/4.0"' \
#                       cmd='addcol err_peak_flux "sqrt(err_peak_flux_076*err_peak_flux_076 + err_peak_flux_084*err_peak_flux_084 + err_peak_flux_092*err_peak_flux_092 + err_peak_flux_099*err_peak_flux_099)"' \
