# Copyright Tobias Wood 2017
# Simple test scripts for QUIT programs

setup() {
load $BATS_TEST_DIRNAME/common.bash
init_tests
# Setup parameters
SPGR_FILE="spgr$EXT"
SPGR_SEQ="\"SPGR\": { \"TR\": 0.01, \"FA\": [3, 8, 13, 18] }"
SIZE="16,16,16"
NOISE="0.001"
TOL="35"
cat > despot1.json <<OUT
{
    $SPGR_SEQ
}
OUT
cat > simulate_despot1.json <<OUT
{
    "PDFile": "PD$EXT",
    "T1File": "T1$EXT",
    $SPGR_SEQ
}
OUT
}

@test "DESPOT1-Simulate" {
qinewimage --size "$SIZE" -g "1 0.8 1.0" PD$EXT
[ -e PD$EXT ]
qinewimage --size "$SIZE" -g "0 0.5 1.5" T1$EXT
[ -e T1$EXT ]
qidespot1 --verbose --simulate=$NOISE $SPGR_FILE < simulate_despot1.json
[ -e $SPGR_FILE ]
}

@test "DESPOT1-LLS" {
qidespot1 --verbose $SPGR_FILE < despot1.json
qidiff --baseline=T1$EXT --input=D1_T1$EXT --noise=$NOISE --tolerance=$TOL --verbose
}

@test "DESPOT1-WLLS" {
qidespot1 --verbose --algo=w $SPGR_FILE < despot1.json
qidiff --baseline=T1$EXT --input=D1_T1$EXT --noise=$NOISE --tolerance=$TOL --verbose
}

@test "DESPOT1-NLLS" {
qidespot1 --verbose --algo=n $SPGR_FILE < despot1.json
qidiff --baseline=T1$EXT --input=D1_T1$EXT --noise=$NOISE --tolerance=$TOL --verbose
}
