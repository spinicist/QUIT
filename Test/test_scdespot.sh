#!/bin/bash

# Tobias Wood 2015
# Simple test script for DESPOT programs

# Tests whether programs run successfully on toy data

# First, create input data

source ./test_common.sh
SILENCE_TESTS="1"

DATADIR="scdespot"
rm -r $DATADIR
mkdir -p $DATADIR
cd $DATADIR

DIMS="11 11 11"
VOXDIMS="2 2 2"
$QUITDIR/qinewimg -d "$DIMS" -v "$VOXDIMS" -f 1 PD.nii
$QUITDIR/qinewimg -d "$DIMS" -v "$VOXDIMS" -g "0 0.5 5" T1.nii
$QUITDIR/qinewimg -d "$DIMS" -v "$VOXDIMS" -g "1 0.05 0.5" T2.nii
$QUITDIR/qinewimg -d "$DIMS" -v "$VOXDIMS" -g "2 -25.0 25.0" f0.nii
$QUITDIR/qinewimg -d "$DIMS" -v "$VOXDIMS" -g "2 0.5 1.5" B1.nii

# Setup parameters
SPGR_FILE="spgr.nii"
SPGR_PAR="3
5 10 15
0.01"
SSFP_FILE="ssfp.nii"
SSFP_FLIP="15 30 45"
SSFP_TR="0.005"
SSFP_PAR="3
$SSFP_FLIP
4
0 90 180 270
$SSFP_TR"
D2GS_PAR="3
$SSFP_FLIP
$SSFP_TR"
MPRAGE_FILE="mprage.nii"
MPRAGE_PAR="5
0.01
11
1
0.5
0"
HIFI_PAR="$SPGR_PAR
$MPRAGE_PAR"

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
END"
echo "$MCSIG_INPUT" > mcsignal.in
run_test "CREATE_REAL_SIGNALS" $QUITDIR/mcsignal --1 -n < mcsignal.in

run_test "CREATE_COMPLEX_SIGNALS" $QUITDIR/mcsignal --1 -x -v <<END_IN
PD.nii
T1.nii
T2.nii
f0.nii
B1.nii
SSFP
$SSFP_PAR
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

run_test "DESPOT1" $QUITDIR/despot1 $SPGR_FILE -n -bB1.nii < despot1.in
compare_test "DESPOT1" T1.nii D1_T1.nii 0.01
run_test "DESPOT1LM" $QUITDIR/despot1 $SPGR_FILE -n -an -oN -bB1.nii < despot1.in
compare_test "DESPOT1LM" T1.nii ND1_T1.nii 0.01
run_test "DESPOT1HIFI" $QUITDIR/despot1hifi $SPGR_FILE $MPRAGE_FILE -M -n < despot1hifi.in
compare_test "HIFI_T1" T1.nii HIFI_T1.nii 0.01
run_test "SSFPGS" $QUITDIR/ssfpgs ssfp_x.nii
run_test "SSFPGSMAG" $QUITDIR/qcomplex -x ssfp_x_lreg.nii -om ssfp_x_lreg_mag.nii
run_test "DESPOT2GS" $QUITDIR/despot2 -e D1_T1.nii ssfp_x_lreg_mag.nii -n -bB1.nii < despot2gs.in
compare_test "DESPOT2GS" T2.nii D2_T2.nii 0.01
run_test "SSFPGS2P" $QUITDIR/ssfpgs -2 ssfp_x.nii
run_test "SSFPGS2PMAG" $QUITDIR/qcomplex -x ssfp_x_lreg_2p.nii -om ssfp_x_lreg_2p_mag.nii
run_test "DESPOT2GS2P" $QUITDIR/despot2 -e D1_T1.nii ssfp_x_lreg_2p_mag.nii -n -bB1.nii -o 2p < despot2gs.in
compare_test "DESPOT2GS2P" T2.nii 2pD2_T2.nii 0.05
run_test "DESPOT2FM" $QUITDIR/despot2fm D1_T1.nii $SSFP_FILE -n -v -S1 -bB1.nii < despot2fm.in
compare_test "DESPOT2FM" T2.nii FM_T2.nii 0.01

cd ..
SILENCE_TESTS="0"
