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
    -m      Specify a mask
    -v      Verbose mode (display progress)
"
MASK=""
VERBOSE=""
while getopts "m:v" opt; do
    case $opt in
        m) MASK="-mas $OPTARG";;
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
MPVOL="dmt0000"
PMVOL="dmt0001"
REFVOL="dmt0002"
PVOL="dmt0003"
MVOL="dmt0004"

log "Calculating ratios"
fslmaths $PVOL -add $MVOL -sub 1 -div 2 -div $REFVOL -mul -100 $MASK DMT_MTR
fslmaths $PMVOL -add $MPVOL -sub 1 -div 2 -div $REFVOL -mul -100 $MASK DMT_eMTR
fslmaths DMT_eMTR -sub DMT_MTR $MASK DMT_dcMTR
fslmaths $PVOL -sub $MVOL -div $REFVOL -mul 100 $MASK DMT_mtasym

log "Cleaning up"
rm $MPVOL.nii* $PMVOL.nii* $REFVOL.nii* $PVOL.nii* $MVOL.nii*
