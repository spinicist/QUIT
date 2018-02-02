#!/bin/bash -eux
##
## qi_composer.sh
##
## Copyright Tobias C Wood 2017 tobias.wood@kcl.ac.uk
##
##  This file is subject to the terms of the Mozilla Public
##  License, v. 2.0. If a copy of the MPL was not distributed with this
##  file, You can obtain one at http://mozilla.org/MPL/2.0/.
##

USAGE="Usage: $0 -r reference_file [opts] input_1 [input_2 ...]

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
    -r      Reference (short echo time) file. REQUIRED
    -t      Keep temporary files
    -m      Output magnitude/phase instead of complex final file
    -v      Verbose mode (display progress)
"

KEEP_TEMP=""
OUTPUT_MAG=""
VERBOSE=""
EXT=".nii.gz"
while getopts "mr:tv" opt; do
    case $opt in
        m) OUTPUT_MAG="1";;
        r) SER="$OPTARG";;
        t) KEEP_TEMP="1";;
        v) VERBOSE="-v";;
    esac
done
shift $(( $OPTIND - 1 ))

if [ -z ${SER-} ]; then
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

while test ${#} -gt 0; do
    IMG="${1}"
    IMG_ROOT="$( basename ${IMG%%.*} )"
    SER_ROOT="$( basename ${SER%%.*} )"
    log "Processing $IMG_ROOT"

    log "Creating magnitude images for registration"
    qicomplex --realimag $IMG -M $TEMP/${IMG_ROOT}_mag${EXT}
    qicomplex --realimag $SER -M $TEMP/${SER_ROOT}_mag${EXT}

    antsMotionCorr -d 3 -a $TEMP/${IMG_ROOT}_mag${EXT} -o $TEMP/${IMG_ROOT}_avg${EXT}
    antsMotionCorr -d 3 -a $TEMP/${SER_ROOT}_mag${EXT} -o $TEMP/${SER_ROOT}_avg${EXT}

    log "Registering reference to input"
    qimask $TEMP/${IMG_ROOT}_avg${EXT} --fillh=2
    ImageMath 3 $TEMP/mask${EXT} MD $TEMP/${IMG_ROOT}_avg_mask${EXT} 3

    antsRegistration --dimensionality 3 --float 0 --interpolation Linear \
        --output [$TEMP/reg,$TEMP/regWarped${EXT}] -x $TEMP/mask${EXT} \
        --initial-moving-transform [$TEMP/${IMG_ROOT}_avg${EXT}, $TEMP/${SER_ROOT}_avg${EXT}, 1] \
        --transform Rigid[0.1] --metric MI[$TEMP/${IMG_ROOT}_avg${EXT}, $TEMP/${SER_ROOT}_avg${EXT}, 1, 32, Regular, 0.25] \
        --convergence [1000x500x250,1e-6,10] --shrink-factors 8x4x2 --smoothing-sigmas 4x2x1vox

    log "Resampling reference"
    antsApplyTransforms --dimensionality 3 --input-image-type 3 \
        --input $SER --reference-image $TEMP/${IMG_ROOT}_avg${EXT} \
        --output $TEMP/${SER_ROOT}_resamp${EXT} \
        --transform $TEMP/reg0GenericAffine.mat --float

    qicomplex --realimag $IMG -X $TEMP/${IMG_ROOT}_x${EXT}
    qicomplex --realimag $TEMP/${SER_ROOT}_resamp${EXT} -P $TEMP/${SER_ROOT}_ph${EXT}

    log "Combining coil images"
    qi_coil_combine --verbose $TEMP/${IMG_ROOT}_x${EXT} --composer=$TEMP/${SER_ROOT}_ph${EXT} --out ${IMG_ROOT}
    if [ -n "$OUTPUT_MAG" ]; then
        log "Writing magnitude/phase output"
        qicomplex -x ${IMG_ROOT}_combined${EXT} -M ${IMG_ROOT}_combined_mag${EXT} -P ${IMG_ROOT}_combined_ph${EXT}
        rm ${IMG_ROOT}_combined${EXT}
    fi
    shift 1
done

if [ -z "$KEEP_TEMP" ]; then
    log "Removing temporary files"
    rm -r $TEMP
fi
