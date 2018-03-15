# Combine all the filtered red catalogues with the priorised fitting results

$java stilts tcat \
in=Week1_red_lownoise_ddmod_rescaled_subtracted_sf_filtered.fits \
in=Week2_red_lownoise_ddmod_rescaled_subtracted_sf_filtered.fits \
in=Week3_red_lownoise_ddmod_rescaled_subtracted_sf_filtered.fits \
in=Week4_red_lownoise_ddmod_rescaled_subtracted_sf_filtered.fits \
out=red_all.fits

for week in Week1 Week2 Week3 Week4
do
#    if [[ -e ${week}_subs.fits ]]; then continue;fi;
    stilts tmatchn multimode=pairs nin=4 matcher=exact \
        in1=${week}_red_lownoise_ddmod_rescaled_subtracted_sf_filtered.fits  values1='uuid'   suffix1='_088' \
        in2=${week}_green_prior_comp.vot values2='uuid' suffix2='_118' \
        in3=${week}_blue_prior_comp.vot values3='uuid' suffix3='_154' \
        in4=${week}_white_prior_comp.vot values4='uuid' suffix4='_200' \
        out=${week}_subs.fits ofmt=fits-basic
done

stilts tcat in=Week1_subs.fits in=Week2_subs.fits in=Week3_subs.fits in=Week4_subs.fits out=Week1-4_subs_long.fits

