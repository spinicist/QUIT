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

USAGE="Usage: $0 input_file reference_file

This script is an implementation of the COMPOSER algorithm (Robinson et al MRM
2017). The actual coil-combination step is implemented in qi_coil_combine. The
remainder of the algorithm involves registration and resampling, which can be
accomplished better using ANTs.

Hence this script requires ANTs.

The inputs are the multi-coil input image and the short echo-time reference
(SER) file. For now, these are assumed to be in Bruker 'complex' format - all
real volumes followed by all complex volumes. The number of coils/volumes in the
two files must match.
"

if [ $# -ne 2 ]; then
echo "$USAGE"
exit 1;
fi

IMG="$1"
SER="$2"

IMG_ROOT="$( basename ${IMG%%.*})"
SER_ROOT="$( basename ${SER%%.*})"

TEMP="composer_working_dir"
mkdir -p $TEMP

EXT="nii.gz"

## First, do magnitude only combination to get images for registration
qicomplex --realimag $IMG -M $TEMP/${IMG_ROOT}_mag.$EXT
qicomplex --realimag $SER -M $TEMP/${SER_ROOT}_mag.$EXT

antsMotionCorr -d 3 -a $TEMP/${IMG_ROOT}_mag.$EXT -o $TEMP/${IMG_ROOT}_avg.$EXT
antsMotionCorr -d 3 -a $TEMP/${SER_ROOT}_mag.$EXT -o $TEMP/${SER_ROOT}_avg.$EXT

## Register SER to input
qimask $TEMP/${IMG_ROOT}_avg.$EXT --fillh=2
ImageMath 3 $TEMP/mask.$EXT MD $TEMP/${IMG_ROOT}_avg_mask.$EXT 3

antsRegistration --verbose 1 --dimensionality 3 --float 0 --interpolation Linear \
    --output [$TEMP/reg,$TEMP/regWarped.$EXT] -x $TEMP/mask.$EXT \
    --initial-moving-transform [$TEMP/${IMG_ROOT}_avg.$EXT, $TEMP/${SER_ROOT}_avg.$EXT, 1] \
    --transform Rigid[0.1] --metric MI[$TEMP/${IMG_ROOT}_avg.$EXT, $TEMP/${SER_ROOT}_avg.$EXT, 1, 32, Regular, 0.25] \
    --convergence [1000x500x250,1e-6,10] --shrink-factors 8x4x2 --smoothing-sigmas 4x2x1vox

## Resample SER to input
antsApplyTransforms --dimensionality 3 --input-image-type 3 \
    --input $SER --reference-image $TEMP/${IMG_ROOT}_avg.$EXT \
    --output $TEMP/${SER_ROOT}_resamp.$EXT \
    --transform $TEMP/reg0GenericAffine.mat --float --verbose

## Now convert to 'real' complex
qicomplex --realimag $IMG -X $TEMP/${IMG_ROOT}_x.$EXT
qicomplex --realimag $TEMP/${SER_ROOT}_resamp.$EXT -X $TEMP/${SER_ROOT}_x.$EXT

## Coil combine
qi_coil_combine $TEMP/${IMG_ROOT}_x.$EXT $TEMP/${SER_ROOT}_x.$EXT --out ${IMG_ROOT}_
qicomplex -x ${IMG_ROOT}_combined.$EXT -M ${IMG_ROOT}_combined_mag.$EXT -P ${IMG_ROOT}_combined_phase.$EXT
