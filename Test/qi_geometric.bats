# Copyright Tobias Wood 2017
# Simple test scripts for QUIT programs

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "Geometric Solution" {

# First, create input data
SIZE="16,16,51"
qinewimage PD.nii --size "$SIZE" -f "1.0"
qinewimage T1.nii --size "$SIZE" -g "0 1.0 2.0"
qinewimage T2.nii --size "$SIZE" -g "1 0.025 0.25"
qinewimage f0.nii --size "$SIZE" -g "2 -200.0 200.0"
qinewimage B1.nii --size "$SIZE" -f 1.0

# Setup parameters
SSFP_FILE="ssfp.nii.gz"
GS_FILE="gs_ref.nii.gz"
SSFP_TR="0.005"
SSFP_FLIP="45"
SSFP_PC="180 270 0 90"
NOISE="2e-3"
qisignal -v -n -x --noise=$NOISE -T1 << END_SIG
PD.nii
T1.nii
T2.nii
f0.nii

$SSFP_FILE
SSFP_ECHO
$SSFP_FLIP
$SSFP_PC
$SSFP_TR
$GS_FILE
SSFP_GS
$SSFP_FLIP
$SSFP_TR
END
END_SIG

qissfpbands ${SSFP_FILE} --ph-incs=4 --verbose
[ -e ssfp_GSL.nii.gz ]
qicomplex -x ssfp_GSL.nii.gz -M gs_mag.nii.gz -P gs_ph.nii.gz --verbose
[ -e gs_mag.nii.gz ]
[ -e gs_ph.nii.gz ]
qicomplex -x gs_ref.nii.gz -M gs_ref_mag.nii.gz --verbose
[ -e gs_ref_mag.nii.gz ]
qidiff --baseline=gs_ref_mag.nii.gz --input=gs_mag.nii.gz --noise=$NOISE --tolerance=50 --verbose

}
