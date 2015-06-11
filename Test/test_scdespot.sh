#!/bin/bash

# Tobias Wood 2015
# Simple test script for DESPOT programs

# Tests whether programs run successfully on toy data

# First, create input data

source ./test_common.sh
SILENCE_TESTS="1"

DATADIR="scdespot"
rm -rf $DATADIR
mkdir -p $DATADIR
cd $DATADIR

DIMS="16 16 25"
VOXDIMS="2 2 2"
$QUITDIR/qinewimage -d "$DIMS" -v "$VOXDIMS" -g "2 1.2 0.8" PD.nii
$QUITDIR/qinewimage -d "$DIMS" -v "$VOXDIMS" -g "0 0.5 2" T1.nii
$QUITDIR/qinewimage -d "$DIMS" -v "$VOXDIMS" -g "1 0.05 0.5" T2.nii
$QUITDIR/qinewimage -d "$DIMS" -v "$VOXDIMS" -g "2 -25.0 25.0" f0.nii
$QUITDIR/qinewimage -d "$DIMS" -v "$VOXDIMS" -g "2 0.75 1.25" B1.nii

# Setup parameters
SPGR_FILE="spgr.nii"
SPGR_PAR="5 10 15 #Test comment
#Test line below has trailing whitespace
0.01 "
SSFP_FILE="ssfp.nii"
SSFP_PAR="15 30 45
0 180
0.005"
SSFP_X_PAR="15 30 45
0 90 180 270
0.005"
D2GS_PAR="15 30 45
0.005"
MPRAGE_FILE="mprage.nii"
MPRAGE_PAR="5
0.01
11
0.5
0"
HIFI_PAR="$SPGR_PAR
$MPRAGE_PAR"
AFI_FILE="afi.nii"
AFI_PAR="55.
0.02 0.1"
# Create input for Single Component
MCSIG_INPUT="PD.nii
T1.nii
T2.nii
f0.nii
B1.nii
SPGR
$SPGR_PAR
$SPGR_FILE
SSFP
$SSFP_PAR
$SSFP_FILE
MPRAGE
$MPRAGE_PAR
$MPRAGE_FILE
AFI
$AFI_PAR
$AFI_FILE
END"
echo "$MCSIG_INPUT" > signal.in
run_test "CREATE_REAL_SIGNALS" $QUITDIR/qisignal --1 -n -v < signal.in

run_test "CREATE_COMPLEX_SIGNALS" $QUITDIR/qisignal --1 -n -x <<END_IN
PD.nii
T1.nii
T2.nii
f0.nii
B1.nii
SSFP
$SSFP_X_PAR
ssfp_x.nii
SSFPEllipse
$D2GS_PAR
ssfp_gs.nii
END
END_IN

echo "$SPGR_PAR" > despot1.in
echo "$HIFI_PAR" > despot1hifi.in
echo "$SSFP_PAR" > despot2fm.in
echo "$D2GS_PAR" > despot2gs.in

run_test "DESPOT1" $QUITDIR/qidespot1 $SPGR_FILE -n -bB1.nii < despot1.in
compare_test "DESPOT1" T1.nii D1_T1.nii 0.01
run_test "DESPOT1LM" $QUITDIR/qidespot1 $SPGR_FILE -n -an -oN -bB1.nii < despot1.in
compare_test "DESPOT1LM" T1.nii ND1_T1.nii 0.01
run_test "DESPOT1HIFI" $QUITDIR/qidespot1hifi $SPGR_FILE $MPRAGE_FILE -M -n -T1 < despot1hifi.in
compare_test "HIFI_T1" T1.nii HIFI_T1.nii 0.01
run_test "AFI" $QUITDIR/qiafi $AFI_FILE
compare_test "AFI_B1" B1.nii AFI_B1.nii 0.01
run_test "SSFPGS" $QUITDIR/qissfpbands ssfp_x.nii
run_test "SSFPGSMAG" $QUITDIR/qicomplex -x ssfp_x_lreg.nii -om ssfp_x_lreg_mag.nii
run_test "DESPOT2GS" $QUITDIR/qidespot2 -e D1_T1.nii ssfp_x_lreg_mag.nii -n --B1 B1.nii < despot2gs.in
compare_test "DESPOT2GS" T2.nii D2_T2.nii 0.01
run_test "SSFPGS2P" $QUITDIR/qissfpbands -2 ssfp_x.nii
run_test "SSFPGS2PMAG" $QUITDIR/qicomplex -x ssfp_x_lreg_2p.nii -om ssfp_x_lreg_2p_mag.nii
run_test "DESPOT2GS2P" $QUITDIR/qidespot2 -e D1_T1.nii ssfp_x_lreg_2p_mag.nii -n -bB1.nii -o 2p < despot2gs.in
compare_test "DESPOT2GS2P" T2.nii 2pD2_T2.nii 0.05

run_test "DESPOT2FMLM" $QUITDIR/qidespot2fm T1.nii $SSFP_FILE -n -S1 -bB1.nii -v -o LM < despot2fm.in
compare_test "DESPOT2FMLM" T2.nii LMFM_T2.nii 0.001
run_test "DESPOT2FMSRC" $QUITDIR/qidespot2fm T1.nii $SSFP_FILE -n -S1 -as -bB1.nii -v -o SRC < despot2fm.in
compare_test "DESPOT2FMSRC" T2.nii SRCFM_T2.nii 0.01

cd ..
SILENCE_TESTS="0"
