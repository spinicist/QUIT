#!/bin/bash -eux

# Tobias Wood 2016
# Simple test scripts for QUIT programs

source ./test_common.sh
SILENCE_TESTS="0"

DATADIR="poly"
mkdir -p $DATADIR
cd $DATADIR
if [ "$(ls -A ./)" ]; then
    rm *
fi

SIZE="32 32 32"
$QUITDIR/qinewimage --size "$SIZE" -f "1" dummy.nii

qipolygen dummy.nii gen.nii --order=2 -v << END_INPUT
0 1 1 1 1 1 1 1 1 1
END_INPUT
qipolyfit gen.nii --order=2 -v

cd ..
SILENCE_TESTS="0"
