#!/bin/bash

if [[ $1 ]] && [[ $2 ]]
then
    GLEAMfits=$1
    if [[ -e $GLEAMfits ]] && [[ -e $2 ]]
    then
        catalogue=$2
        GLEAM=`echo $GLEAMfits | sed "s/.fits//"`

        zero.py ${GLEAMfits} ${GLEAM}_zero.fits

        # Select and restore the bit of the EGC catalogue that overlaps
        ./select_subset.py --catalogue=$catalogue --output=${GLEAM}_subset.fits --postage=$GLEAMfits # Why is this script launched differently?
        AeRes -c ${GLEAM}_subset.fits -f ${GLEAM}_zero.fits -r ${GLEAM}_restore.fits --add --mask --sigma=0.5
        rm ${GLEAM}_zero.fits 

        gbmaj=`pyhead.py -p BMAJ $GLEAMfits | awk '{print 3600*$3}'`
        gbmin=`pyhead.py -p BMIN $GLEAMfits | awk '{print 3600*$3}'`
        gbpa=`pyhead.py -p BPA $GLEAMfits | awk '{print $3}'`

        subtract.py $GLEAMfits ${GLEAM}_restore.fits ${GLEAM}_subtracted.fits
    else
        echo "$1 or $2 doesn't exist."
    fi
else
    echo "need to specify fits file and catalogue"
    exit 1
fi

exit 0
