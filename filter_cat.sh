#!/bin/bash

# Filter catalogues to avoid including a bunch of unnecessary data in the priorised fit

infile=$1
outfile=`echo $infile | sed "s/_comp.fits/_filtered.fits/"`

stilts tpipe in=$infile cmd='addskycoords -inunit deg fk5 gal ra dec GAL_LONG GAL_LAT' \
                       cmd='select "((dec<=30)&&(abs(GAL_LAT)>11)&&(err_ra!=-1.&&err_a!=-1.))"' \
                       cmd='select "!((dec>=0)&&(ra>8.*22.))"' \
                       cmd='select "skyDistanceDegrees(ra,dec,hmsToDegrees(05,23,34.6),dmsToDegrees(-69,45,22)) > 6"' \
                       cmd='select "skyDistanceDegrees(ra,dec,hmsToDegrees(00,52,38.0),dmsToDegrees(-72,48,1)) > 3"' \
                       cmd='select "skyDistanceDegrees(ra,dec,hmsToDegrees(05,19,49.7),dmsToDegrees(-45,46,53.9)) > 1"' \
                       cmd='select "skyDistanceDegrees(ra,dec,hmsToDegrees(16,51,11.4),dmsToDegrees(4,59,20)) > 1"' \
                       cmd='select "skyDistanceDegrees(ra,dec,hmsToDegrees(12,30,49.4),dmsToDegrees(12,23,28)) > 1"' \
                       cmd='select "skyDistanceDegrees(ra,dec,hmsToDegrees(13,25,27.6),dmsToDegrees(-43,1,9)) > 10"' \
                       cmd='keepcols "ra err_ra dec err_dec int_flux err_int_flux peak_flux err_peak_flux a err_a b err_b pa err_pa local_rms uuid"' \
                       out=temp.fits

if [[ ${1:0:5} == "Week1" ]]
then
    stilts tpipe in=temp.fits cmd='select "(ra>=21*15.)"'\
    out=$outfile
    rm temp.fits
elif [[ ${1:0:5} == "Week2" ]]
then
    stilts tpipe in=temp.fits cmd='select "((ra>=0)&&(ra<8.*15.))"'\
    out=$outfile
    rm temp.fits
elif [[ ${1:0:5} == "Week3" ]]
then
    stilts tpipe in=temp.fits cmd='select "((ra>=8*15.)&&(ra<14.5*15.))"'\
    out=$outfile
    rm temp.fits
elif [[ ${1:0:5} == "Week4" ]]
then
    stilts tpipe in=temp.fits cmd='select "((ra>=14.5*15.)&&(ra<21*15.))"'\
    out=$outfile
    rm temp.fits
else
    mv $outfile temp.fits
fi


