#!/bin/bash -e

# Tobias Wood 2015
# Test script for ssfpbands etc.

source ./test_common.sh
SILENCE_TESTS="1"

# First, create input data
DATADIR="ssfp"
mkdir -p $DATADIR
cd $DATADIR
if [ "$(ls -A ./)" ]; then
    rm *
fi

DIMS="16 16 25"
VOXDIMS="2 2 2"
$QUITDIR/qinewimage -d "$DIMS" -v "$VOXDIMS" -f "1.0" PD.nii
$QUITDIR/qinewimage -d "$DIMS" -v "$VOXDIMS" -f "1.0" T1.nii
$QUITDIR/qinewimage -d "$DIMS" -v "$VOXDIMS" -f "0.1" T2.nii
$QUITDIR/qinewimage -d "$DIMS" -v "$VOXDIMS" -g "2 -25.0 25.0" f0.nii
$QUITDIR/qinewimage -d "$DIMS" -v "$VOXDIMS" -g "2 0.75 1.25" B1.nii

SSFP_FILE="ssfp.nii"
SSFP_OUT_FILE="ssfp_lreg.nii"
SSFP_MAG_FILE="ssfp_lreg_mag.nii"
SSFP_PAR="15 30 45
0 90 180 270
0.005"
GS_FILE="gs.nii"
GS_MAG_FILE="gs_mag.nii"
GS_PAR="15 30 45
0.005"

run_test "CREATE_COMPLEX_SIGNALS" $QUITDIR/qisignal --1 -n -x <<END_IN
PD.nii
T1.nii
T2.nii
f0.nii
B1.nii
$SSFP_FILE
SSFP
$SSFP_PAR
$GS_FILE
SSFPEllipse
$GS_PAR
END
END_IN

run_test "SSFPGS" $QUITDIR/qissfpbands $SSFP_FILE
$QUITDIR/qicomplex -x $SSFP_OUT_FILE -M $SSFP_MAG_FILE
$QUITDIR/qicomplex -x $GS_FILE -M $GS_MAG_FILE
compare_test "DESPOT2GS" $GS_MAG_FILE $SSFP_MAG_FILE 0.01
