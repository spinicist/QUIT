#!/bin/bash -eux

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

SIZE="25 25 51"
$QUITDIR/qinewimage --size "$SIZE" -f "1.0" PD.nii
$QUITDIR/qinewimage --size "$SIZE" -f "1.0" T1.nii
$QUITDIR/qinewimage --size "$SIZE" -f "0.1" T2.nii #-g "0 0.02 0.25" T2.nii
$QUITDIR/qinewimage --size "$SIZE" -g "2 -300.0 300.0" f0.nii
$QUITDIR/qinewimage --size "$SIZE" -f "1.0" B1.nii #-g "2 0.5 1.0"

SSFP_FILE="ssfp.nii"
SSFP_FLIP="20 60"
SSFP_PINC="180 0 90 270"
SSFP_TR="0.005"
GS_REF_FILE="gs_ref.nii"
GS_REF_MAG="gs_ref_mag.nii"
GS_FILE="ssfp_GS.nii"
GSM_FILE="ssfp_GSM.nii"
GSL_FILE="ssfp_GSL.nii"
GS2_FILE="ssfp_GSL2.nii"
CS_FILE="ssfp_CS.nii"
MAG_FILE="ssfp_MagMean.nii"
RMS_FILE="ssfp_RMS.nii"
MAX_FILE="ssfp_Max.nii"
NOISE="0.001"
run_test "CREATE_COMPLEX_SIGNALS" $QUITDIR/qisignal --model=1 -n --noise=$NOISE -x <<END_IN
PD.nii
T1.nii
T2.nii
f0.nii
B1.nii
$SSFP_FILE
SSFP_ECHO
$SSFP_FLIP
$SSFP_PINC
$SSFP_TR
$GS_REF_FILE
SSFP_GS
$SSFP_FLIP
$SSFP_TR
END
END_IN
$QUITDIR/qicomplex -x $GS_REF_FILE -M $GS_REF_MAG -P gs_ref_ph.nii
fslmaths gs_ref_ph -div 3.141592 -div $SSFP_TR gs_ref_f0.nii

$QUITDIR/qicomplex -x $SSFP_FILE -M ssfp_mag.nii -P ssfp_ph.nii
run_test "GS"  $QUITDIR/qissfpbands $SSFP_FILE --regularise N --alt-order
$QUITDIR/qicomplex -x $GS_FILE -M gs_mag.nii -P gs_ph.nii
fslmaths gs_ph -div 3.141592 -div $SSFP_TR gs_f0.nii

$QUITDIR/qidespot2fm -n -v -bB1.nii T1.nii ssfp_mag.nii --asym <<END_IN
$SSFP_FLIP
$SSFP_PINC
$SSFP_TR
END_IN

#run_test "GSM" $QUITDIR/qissfpbands $SSFP_FILE --regularise M      --magnitude -o $GSM_FILE
#run_test "GSL" $QUITDIR/qissfpbands $SSFP_FILE --regularise L      --magnitude -o $GSL_FILE
#run_test "GS2" $QUITDIR/qissfpbands $SSFP_FILE --regularise L -2   --magnitude -o $GS2_FILE
#run_test "CS"  $QUITDIR/qissfpbands $SSFP_FILE --cs      --magnitude -o $CS_FILE
#run_test "MAG" $QUITDIR/qissfpbands $SSFP_FILE --magmean --magnitude -o $MAG_FILE
#run_test "RMS" $QUITDIR/qissfpbands $SSFP_FILE ---regularisems     --magnitude -o $RMS_FILE
#run_test "MAX" $QUITDIR/qissfpbands $SSFP_FILE --max     --magnitude -o $MAX_FILE

#compare_test "GS"  $GS_FILE  $GS_REF_MAG $NOISE 50
#compare_test "GSM" $GSM_FILE $GS_REF_MAG $NOISE 50
#compare_test "GSL" $GSL_FILE $GS_REF_MAG $NOISE 50
#compare_test "GS2" $GS2_FILE $GS_REF_MAG $NOISE 50
#compare_test "CS"  $CS_FILE  $GS_REF_MAG $NOISE 50
#compare_test "MAG" $MAG_FILE $GS_REF_MAG $NOISE 50
#compare_test "RMS" $RMS_FILE $GS_REF_MAG $NOISE 50
#compare_test "MAX" $MAX_FILE $GS_REF_MAG $NOISE 50
