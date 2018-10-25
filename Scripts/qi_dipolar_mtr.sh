#!/bin/bash -eu
##
## qi_dipolar_mtr.sh
##
## Copyright Tobias C Wood 2018 tobias.wood@kcl.ac.uk
##
##  This file is subject to the terms of the Mozilla Public
##  License, v. 2.0. If a copy of the MPL was not distributed with this
##  file, You can obtain one at http://mozilla.org/MPL/2.0/.
##

USAGE="Usage: $0 [opts] input

This script calculates MTR / dcMTR from a time-series with 5 volumes:
1 : +/- saturation
2 : -/+ saturation
3 : No saturation
4 : + saturation
5 : - saturation

Use FSL

Options:
    -m MASK  Specify a mask
    -o ORDER Change order of images REF P M PM MP
    -v       Verbose mode (display progress)
"
MASK=""
ORDER=(2 3 4 0 1)
VERBOSE=""
while getopts "m:o:v" opt; do
    case $opt in
        m) MASK="-mas $OPTARG";;
        o) ORDER=($OPTARG);;
        v) VERBOSE="-v";;
    esac
done
shift $(( $OPTIND - 1 ))

if [ $# -ne 1 ]; then
    echo "$USAGE"
    exit 1;
fi

function log() {
    if [ -n "$VERBOSE" ]; then
        echo "$1"
    fi
}

INPUT=$1

log "Splitting input"
fslsplit $INPUT dmt -t
SPLIT_IMGS=("dmt0000" "dmt0001" "dmt0002" "dmt0003" "dmt0004")
REFVOL=${SPLIT_IMGS[${ORDER[0]}]}
PVOL=${SPLIT_IMGS[${ORDER[1]}]}
MVOL=${SPLIT_IMGS[${ORDER[2]}]}
MPVOL=${SPLIT_IMGS[${ORDER[3]}]}
PMVOL=${SPLIT_IMGS[${ORDER[4]}]}

log "Calculating ratios"
fslmaths $PVOL -add $MVOL -div 2 -div $REFVOL -sub 1 -mul -100 -thr -100 -uthr 100 $MASK DMT_MTR
fslmaths $PMVOL -add $MPVOL -div 2 -div $REFVOL -sub 1 -mul -100 -thr -100 -uthr 100 $MASK DMT_eMTR
fslmaths DMT_eMTR -sub DMT_MTR -thr -100 -uthr 100 $MASK DMT_dcMTR
fslmaths $PVOL -sub $MVOL -div $REFVOL -mul 100 -thr -100 -uthr 100 $MASK DMT_mtasym

log "Cleaning up"
rm $MPVOL.nii* $PMVOL.nii* $REFVOL.nii* $PVOL.nii* $MVOL.nii*
