#!/bin/bash

week=Week2
weeks="Week1 Week2 Week3 Week4"

ncpus=7

for week in $weeks
do
    colors="green blue white"

    redbeam=`pyhead.py -p BMAJ ${week}_red_lownoise_ddmod_rescaled.fits | awk '{print $3}'` # what is pyead.py?

    for color in $colors
    do
        beam=`pyhead.py -p BMAJ ${week}_${color}_lownoise_ddmod_rescaled.fits | awk '{print $3}'`
        psfratio=`echo "$beam / $redbeam" | bc -l`
        echo $psfratio
        aegean --cores=${ncpus} --telescope=mwa --maxsummits=5 --noise=${week}_${color}_lownoise_ddmod_rms_rescaled.fits --background=${week}_${color}_lownoise_ddmod_bkg_rescaled.fits --out=/dev/null --table=${week}_${color}_prior.vot,${week}_${color}_prior.reg --input=${week}_red_lownoise_ddmod_rescaled_subtracted_sf_filtered.fits --priorized=1 --ratio=$psfratio --psf=${week}_${color}_lownoise_comp_psf.fits ${week}_${color}_lownoise_ddmod_rescaled.fits
    done
done
