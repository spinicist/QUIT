#!/bin/bash -eu

# Tobias Wood 2016
# Simple test scripts for QUIT programs

source ./test_common.sh
SILENCE_TESTS="1"

DATADIR="poly"
mkdir -p $DATADIR
cd $DATADIR
if [ "$(ls -A ./)" ]; then
    rm *
fi

SIZE="32 32 32"
$QUITDIR/qinewimage --size "$SIZE" --step="2 0 1 2" mask.nii --verbose

$QUITDIR/qipolyimg mask.nii baseline.nii --order=2 -v << END_INPUT
0 0 0
0 1 2 4 1 0 0 1 0 1
END_INPUT
$QUITDIR/qipolyfit baseline.nii --order=2 --print-terms > unmasked_terms.txt
$QUITDIR/qipolyfit baseline.nii --order=2 --print-terms --mask=mask.nii > masked_terms.txt
$QUITDIR/qipolyimg mask.nii unmasked.nii --order=2 < unmasked_terms.txt
$QUITDIR/qipolyimg mask.nii masked.nii --order=2 < masked_terms.txt
compare_test "UNMASKED" baseline.nii unmasked.nii 1 1
compare_test "MASKED"   baseline.nii masked.nii   1 1
cd ..
SILENCE_TESTS="0"
