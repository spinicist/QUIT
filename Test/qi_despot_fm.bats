# Copyright Tobias Wood 2017
# Tests for DESPOT2-FM

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "DESPOT2-FM" {

SIZE="16,16,17"
qinewimage PD.nii --size "$SIZE" -f 1.0
qinewimage T1.nii --size "$SIZE" -f 1.0
qinewimage T2.nii --size "$SIZE" -g "1 0.025 0.125"
qinewimage f0.nii --size "$SIZE" -g "2 -150.0 150.0"
qinewimage B1.nii --size "$SIZE" -f 1.0

# Setup parameters
SSFP_FILE="ssfp.nii"
SSFP_FLIP="12 65"
SSFP_PINC="0 90 180 270"
SSFP_TR="0.005"
SSFP_Trf="0.001"

NOISE="0.002"
qisignal --model=1 -n -v --noise=$NOISE << END_SIG
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
END
END_SIG

echo "$SSFP_FLIP
$SSFP_PINC
$SSFP_TR" > fm_in.txt

qidespot2fm -n -v -bB1.nii T1.nii ${SSFP_FILE} --asym < fm_in.txt
img_diff T2.nii FM_T2.nii $NOISE 50

}
