# Copyright Tobias Wood 2017
# Simple test scripts for QUIT programs

setup() {
load $BATS_TEST_DIRNAME/common.bash
init_tests
# Setup parameters
SPGR_FILE="spgr$EXT"
MPRAGE_FILE="mprage$EXT"
SPGR=' "SPGR": { "TR": 0.01, "FA": [3, 20] }'
MPRAGE=' "MPRAGE": { "FA": 5, "TR": 0.01, "TI": 0.45, "TD": 0, "eta": 1, "ETL": 64, "k0": 0 }'
SIZE="16,16,16"
NOISE="0.001"
TOL="70"
cat > hifi.json <<OUT
{
    "PDFile": "PD$EXT",
    "T1File": "T1$EXT",
    "B1File": "B1$EXT",
    $SPGR,
    $MPRAGE
}
OUT
}

@test "DESPOT1HIFI-Simulate" {
qinewimage --size "$SIZE" -g "0 0.8 1.0" PD$EXT
[ -e PD$EXT ]
qinewimage --size "$SIZE" -g "1 0.5 1.5" T1$EXT
[ -e T1$EXT ]
qinewimage --size "$SIZE" -g "2 0.5 1.5" B1$EXT
[ -e B1$EXT ]
qidespot1hifi --verbose --simulate=$NOISE $SPGR_FILE $MPRAGE_FILE < hifi.json
}

@test "DESPOT1HIFI" {
qidespot1hifi --verbose $SPGR_FILE $MPRAGE_FILE < hifi.json
qidiff --baseline=T1.nii --input=HIFI_T1.nii --noise=$NOISE --tolerance=$TOL --verbose
qidiff --baseline=B1.nii --input=HIFI_B1.nii --noise=$NOISE --tolerance=$TOL --verbose
}
