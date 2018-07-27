# Copyright Tobias Wood 2017
# Simple test scripts for QUIT programs

setup() {
load $BATS_TEST_DIRNAME/common.bash
init_tests
SSFP_FILE="ssfp$EXT"
SSFP_SEQ="\"SSFP\": { \"TR\": 0.01, \"FA\": [15, 30, 45, 60], \"PhaseInc\": [180, 180, 180, 180] }"
SIZE="16,16,16"
NOISE="0.001"
TOL="20"
cat > despot2.json <<OUT
{
    $SSFP_SEQ
}
OUT
cat > simulate_despot2.json <<OUT
{
    "PDFile": "PD$EXT",
    "T2File": "T2$EXT",
    $SSFP_SEQ
}
OUT
}

@test "DESPOT2-Simulate" {
qinewimage --size "$SIZE" --fill "1" PD$EXT
[ -e PD$EXT ]
qinewimage --size "$SIZE" -g "0 0.8 1.2" T1$EXT
[ -e T1$EXT ]
qinewimage --size "$SIZE" -g "1 0.06 0.1" T2$EXT
[ -e T2$EXT ]
qidespot2 --verbose --simulate=$NOISE T1$EXT $SSFP_FILE < simulate_despot2.json
}

@test "DESPOT2-LLS" {
qidespot2 --verbose T1$EXT $SSFP_FILE < despot2.json
qidiff --baseline=T2$EXT --input=D2_T2$EXT --noise=$NOISE --tolerance=$TOL --verbose
}

@test "DESPOT2-WLLS" {
qidespot2 --verbose --algo=w T1$EXT $SSFP_FILE < despot2.json
qidiff --baseline=T2$EXT --input=D2_T2$EXT --noise=$NOISE --tolerance=$TOL --verbose
}

@test "DESPOT2-NLLS" {
qidespot2 --verbose --algo=n T1$EXT $SSFP_FILE < despot2.json
qidiff --baseline=T2$EXT --input=D2_T2$EXT --noise=$NOISE --tolerance=$TOL --verbose
}

@test "DESPOT2GS-Simulate" {
qinewimage --size "$SIZE" -g "1 0.9 1.1" PD$EXT
[ -e PD$EXT ]
qinewimage --size "$SIZE" -g "0 0.8 1.2" T1$EXT
[ -e T1$EXT ]
qinewimage --size "$SIZE" -g "2 0.05 0.15" T2$EXT
[ -e T2$EXT ]
qidespot2 --verbose --gs --simulate=$NOISE T1$EXT $SSFP_FILE < simulate_despot2.json
}

@test "DESPOT2GS-LLS" {
qidespot2 --verbose --gs T1$EXT $SSFP_FILE < despot2.json
qidiff --baseline=T2$EXT --input=D2_T2$EXT --noise=$NOISE --tolerance=$TOL --verbose
}

@test "DESPOT2GS-WLLS" {
qidespot2 --verbose --gs --algo=w T1$EXT $SSFP_FILE < despot2.json
qidiff --baseline=T2$EXT --input=D2_T2$EXT --noise=$NOISE --tolerance=$TOL --verbose
}

@test "DESPOT2GS-NLLS" {
qidespot2 --verbose --gs --algo=n T1$EXT $SSFP_FILE < despot2.json
qidiff --baseline=T2$EXT --input=D2_T2$EXT --noise=$NOISE --tolerance=$TOL --verbose
}