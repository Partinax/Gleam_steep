#!/bin/bash

weeks="Week1 Week3 Week4"

for week in $weeks
do
    infile=${week}_red_lownoise_ddmod_rescaled_subtracted.fits
    psf=${week}_red_lownoise_comp_psf.fits
    rms=${week}_red_lownoise_ddmod_rms_rescaled.fits
    bkg=${week}_red_lownoise_ddmod_bkg_rescaled.fits

    outfile=`echo $infile | sed "s/.fits/_sf.fits/"`
    outreg=`echo $infile | sed "s/.fits/_sf.reg/"`

    if [[ ! -e no_gp.mim ]]
    then
        MIMAS -g +c 0 -90 80 +c 0 90 80 -o no_gp.mim
    fi

    aegean --cores=7 --out=/dev/null --table=$outfile,$outreg --noise=$rms --psf=$psf --background=$bkg --region=no_gp.mim $infile
done
