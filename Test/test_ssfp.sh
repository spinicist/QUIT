#!/bin/bash -e

# Tobias Wood 2015
# Test script for ssfpbands etc.

source ./test_common.sh
SILENCE_TESTS="0"

# First, create input data
DATADIR="ssfp"
mkdir -p $DATADIR
cd $DATADIR
if [ "$(ls -A ./)" ]; then
    rm *
fi

DIMS="32 32 101"
VOXDIMS="2 2 2"
$QUITDIR/qinewimage -d "$DIMS" -v "$VOXDIMS" -f "1.0" PD.nii
$QUITDIR/qinewimage -d "$DIMS" -v "$VOXDIMS" -f "1.0" T1.nii
$QUITDIR/qinewimage -d "$DIMS" -v "$VOXDIMS" -f "0.1" T2.nii
$QUITDIR/qinewimage -d "$DIMS" -v "$VOXDIMS" -g "2 -100.0 100.0" f0.nii
$QUITDIR/qinewimage -d "$DIMS" -v "$VOXDIMS" -f "1.0" B1.nii

SSFP_FILE="ssfp.nii"
SSFP_MAG="ssfp_mag.nii"
SSFP_PAR="5 10 15 20 25 30 35 40 45 50 55 60
0 90 180 270
0.005"
GS_FILE="gs.nii"
CS_FILE="cs.nii"
SUM_FILE="sum.nii"
SOS_FILE="sos.nii"
MAX_FILE="max.nii"
GS_PAR="5 10 15 20 25 30 35 40 45 50 55 60
0.005"

run_test "CREATE_COMPLEX_SIGNALS" $QUITDIR/qisignal --1 -n -N 0.001 -x <<END_IN
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

run_test "GS" $QUITDIR/qissfpbands $SSFP_FILE -R L -2 --magnitude -o $GS_FILE
run_test "CS" $QUITDIR/qissfpbands $SSFP_FILE --cs --magnitude -o $CS_FILE
run_test "SUM" $QUITDIR/qissfpbands $SSFP_FILE --magsum --magnitude -o $SUM_FILE
run_test "SOS" $QUITDIR/qissfpbands $SSFP_FILE --sos --magnitude -o $SOS_FILE
run_test "MAX" $QUITDIR/qissfpbands $SSFP_FILE --max --magnitude -o $MAX_FILE

$QUITDIR/qicomplex -x $SSFP_FILE -M $SSFP_MAG
compare_test "GS" $GS_FILE $SSFP_MAG 0.01
compare_test "CS" $CS_FILE $SSFP_MAG 0.01
compare_test "SUM" $SUM_FILE $SSFP_MAG 0.01
compare_test "SOS" $SOS_FILE $SSFP_MAG 0.01
compare_test "MAX" $MAX_FILE $SSFP_MAG 0.01
