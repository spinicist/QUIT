#!/bin/bash -e

# Tobias Wood 2015
# Simple test script for FM, it runs slowly so test
# separately to other single-component programs.

# Tests whether programs run successfully on toy data

# First, create input data

source ./test_common.sh
SILENCE_TESTS="1"

DATADIR="fm"
rm -rf $DATADIR
mkdir -p $DATADIR
cd $DATADIR

DIMS="16 10 25"
VOXDIMS="2 2 2"
$QUITDIR/qinewimage PD.nii -d "$DIMS" -v "$VOXDIMS" -g "2 1.2 0.8"
$QUITDIR/qinewimage T1.nii -d "$DIMS" -v "$VOXDIMS" -g "0 0.5 2"
$QUITDIR/qinewimage T2.nii -d "$DIMS" -v "$VOXDIMS" -g "1 0.05 0.5"
$QUITDIR/qinewimage f0.nii -d "$DIMS" -v "$VOXDIMS" -g "2 -50.0 50.0"
$QUITDIR/qinewimage B1.nii -d "$DIMS" -v "$VOXDIMS" -f 1.0

# Setup parameters
SPGR_FILE="spgr.nii"
SPGR_FLIP="5 10 15 #Test comment"
#Test line below has trailing whitespace
SPGR_TR="0.01 "
SPGR_Trf="0.002"
SPGR_TE="0.004"
SSFP_FILE="ssfp.nii"
SSFP_FLIP="5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90"
SSFP_PC="0 180"
SSFP_TR="0.005"
SSFP_Trf="0.0025"

run_test "CREATE_SIGNALS" $QUITDIR/qisignal --1 -n -v -x << END_SIG
PD.nii
T1.nii
T2.nii
f0.nii
B1.nii
SPGR
$SPGR_FLIP
$SPGR_TR
x$SPGR_FILE
SSFP
$SSFP_FLIP
$SSFP_PC
$SSFP_TR
x$SSFP_FILE
SPGRFinite
$SPGR_FLIP
$SPGR_TR
$SPGR_Trf
$SPGR_TE
xF${SPGR_FILE}
SSFPFinite
$SSFP_FLIP
$SSFP_PC
$SSFP_TR
$SSFP_Trf
xF${SSFP_FILE}
END
END_SIG

$QUITDIR/qicomplex -x x$SPGR_FILE -M $SPGR_FILE
$QUITDIR/qicomplex -x x$SSFP_FILE -M $SSFP_FILE -P p$SSFP_FILE
$QUITDIR/qicomplex -x xF$SPGR_FILE -M F$SPGR_FILE
$QUITDIR/qicomplex -x xF$SSFP_FILE -M F$SSFP_FILE -P F$SSFP_FILE

run_test "FM_not_F" $QUITDIR/qidespot2fm T1.nii $SSFP_FILE -n -bB1.nii -v << END_FM
$SSFP_FLIP
$SSFP_PC
$SSFP_TR
END_FM

#run_test "FM_F" $QUITDIR/qidespot2fm T1.nii F${SSFP_FILE} -n -bB1.nii -v -o F --finite << END_FM
#$SSFP_FLIP
#$SSFP_PC
#$SSFP_TR
#$SSFP_Trf
#END_FM

#run_test "FM_F_as_not_F" $QUITDIR/qidespot2fm T1.nii F${SSFP_FILE} -n -bB1.nii -v -o FnF << END_FM
#$SSFP_FLIP
#$SSFP_PC
#$SSFP_TR
#END_FM

compare_test "FM_not_F"      T2.nii FM_T2.nii 0.001
#compare_test "FM_F"          T2.nii FFM_T2.nii 0.001
#compare_test "FM_F_as_not_F" T2.nii FnFFM_T2.nii 0.001

cd ..
SILENCE_TESTS="0"
