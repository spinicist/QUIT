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
SSFP_PAR="5 10 15 20 25 30 35 40 45 50 55 60
0 90 180 270
0.005"
GS_REF_FILE="gs_ref.nii"
GS_REF_MAG="gs_ref_mag.nii"
GS_FILE="gs.nii"
GSM_FILE="gsm.nii"
GSL_FILE="gsl.nii"
GS2_FILE="gs2.nii"
CS_FILE="cs.nii"
MAG_FILE="mag.nii"
RMS_FILE="rms.nii"
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
$GS_REF_FILE
SSFPEllipse
$GS_PAR
END
END_IN
$QUITDIR/qicomplex -x $GS_REF_FILE -M $GS_REF_MAG

run_test "GS"  $QUITDIR/qissfpbands $SSFP_FILE -R N      --magnitude -o $GS_FILE
run_test "GSM" $QUITDIR/qissfpbands $SSFP_FILE -R M      --magnitude -o $GSM_FILE
run_test "GSL" $QUITDIR/qissfpbands $SSFP_FILE -R L      --magnitude -o $GSL_FILE
run_test "GS2" $QUITDIR/qissfpbands $SSFP_FILE -R L -2   --magnitude -o $GS2_FILE
run_test "CS"  $QUITDIR/qissfpbands $SSFP_FILE --cs      --magnitude -o $CS_FILE
run_test "MAG" $QUITDIR/qissfpbands $SSFP_FILE --magmean --magnitude -o $MAG_FILE
run_test "RMS" $QUITDIR/qissfpbands $SSFP_FILE --rms     --magnitude -o $RMS_FILE
run_test "MAX" $QUITDIR/qissfpbands $SSFP_FILE --max     --magnitude -o $MAX_FILE

compare_test "GS"  $GS_FILE  $GS_REF_MAG 0.01
compare_test "GSM" $GSM_FILE $GS_REF_MAG 0.01
compare_test "GSL" $GSL_FILE $GS_REF_MAG 0.01
compare_test "GS2" $GS2_FILE $GS_REF_MAG 0.01
compare_test "CS"  $CS_FILE  $GS_REF_MAG 0.01
compare_test "MAG" $MAG_FILE $GS_REF_MAG 0.01
compare_test "RMS" $RMS_FILE $GS_REF_MAG 0.01
compare_test "MAX" $MAX_FILE $GS_REF_MAG 0.01
