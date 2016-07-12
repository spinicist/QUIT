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
$QUITDIR/qinewimage --size "$SIZE" -f "1" dummy.nii

qipolyimg dummy.nii gen.nii --order=2 -v << END_INPUT
0 1 1 1 1 1 1 1 1 1
END_INPUT
qipolyfit gen.nii --order=2 --print-terms

cd ..
SILENCE_TESTS="0"
