#!/bin/bash -eu
##
## qi_composer.sh
##
## Copyright Tobias C Wood 2017 tobias.wood@kcl.ac.uk
##
##  This file is subject to the terms of the Mozilla Public
##  License, v. 2.0. If a copy of the MPL was not distributed with this
##  file, You can obtain one at http://mozilla.org/MPL/2.0/.
##

USAGE="Usage: $0 [opts] input reference

This script is an implementation of the COMPOSER algorithm (Robinson et al MRM
2017). The actual coil-combination step is implemented in qi_coil_combine. The
remainder of the algorithm involves registration and resampling, which can be
accomplished better using ANTs.

Hence this script requires ANTs.

The inputs are the multi-coil input image and the short echo-time reference
(SER) file. For now, these are assumed to be in Bruker 'complex' format - all
real volumes followed by all complex volumes. The number of coils/volumes in the
two files must match.

Options:
    -m      Output magnitude/phase instead of complex final file
    -o      Specify output prefix (otherwise use input basename)
    -t      Keep temporary files
    -v      Verbose mode (display progress)
"

KEEP_TEMP=""
OUTPUT_MAG=""
VERBOSE=""
EXT=".nii.gz"
IMG_ROOT=""
while getopts "c:mo:tv" opt; do
    case $opt in
        m) OUTPUT_MAG="1";;
        o) IMG_ROOT="$OPTARG";;
        t) KEEP_TEMP="1";;
        v) VERBOSE="-v";;
    esac
done
shift $(( $OPTIND - 1 ))

if [ $# -ne 2 ]; then
    echo "$USAGE"
    exit 1;
fi

function log() {
    if [ -n "$VERBOSE" ]; then
        echo "$1"
    fi
}

log "Creating COMPOSER temporary directory"
TEMP="composer_working_dir"
mkdir -p $TEMP

IMG="$1"
if [ -z ${IMG_ROOT} ]; then
    IMG_ROOT="$( basename ${IMG%%.*} )"
fi
SER="$2"
SER_ROOT="$( basename ${SER%%.*} )"
log "Processing $IMG"

log "Creating magnitude images for registration"
qicomplex --realimag $IMG -M $TEMP/${IMG_ROOT}_mag${EXT}
qicomplex --realimag $SER -M $TEMP/${SER_ROOT}_mag${EXT}

antsMotionCorr -d 3 -a $TEMP/${IMG_ROOT}_mag${EXT} -o $TEMP/${IMG_ROOT}_avg${EXT}
antsMotionCorr -d 3 -a $TEMP/${SER_ROOT}_mag${EXT} -o $TEMP/${SER_ROOT}_avg${EXT}

log "Registering reference to input"
antsRegistration --dimensionality 3 --float 0 --interpolation Linear \
    --output [$TEMP/reg,$TEMP/regWarped${EXT}] \
    --initial-moving-transform [$TEMP/${IMG_ROOT}_avg${EXT}, $TEMP/${SER_ROOT}_avg${EXT}, 1] \
    --transform Rigid[0.1] --metric MI[$TEMP/${IMG_ROOT}_avg${EXT}, $TEMP/${SER_ROOT}_avg${EXT}, 1, 32, Regular, 0.25] \
    --convergence [1000x500x250,1e-6,10] --shrink-factors 8x4x2 --smoothing-sigmas 4x2x1vox

log "Resampling reference"
antsApplyTransforms --dimensionality 3 --input-image-type 3 \
    --input $SER --reference-image $TEMP/${IMG_ROOT}_avg${EXT} \
    --output $TEMP/${SER_ROOT}_resamp${EXT} \
    --transform $TEMP/reg0GenericAffine.mat --float

qicomplex --realimag $IMG -X $TEMP/${IMG_ROOT}_x${EXT}
qicomplex --realimag $TEMP/${SER_ROOT}_resamp${EXT} -X $TEMP/${SER_ROOT}_x${EXT}

log "Combining coil images"
qi_coil_combine $TEMP/${IMG_ROOT}_x${EXT} --verbose --composer=$TEMP/${SER_ROOT}_x${EXT} --out ${IMG_ROOT}
if [ -n "$OUTPUT_MAG" ]; then
    log "Writing magnitude/phase output"
    qicomplex -x ${IMG_ROOT}_combined${EXT} -M ${IMG_ROOT}_combined_mag${EXT} -P ${IMG_ROOT}_combined_ph${EXT}
    if [ -z "$KEEP_TEMP" ]; then
        rm ${IMG_ROOT}_combined${EXT}
    fi
fi

if [ -z "$KEEP_TEMP" ]; then
    log "Removing temporary files"
    rm -r $TEMP
fi
