# Copyright Tobias Wood 2017
# Tests for DESPOT2-FM

setup() {
load $BATS_TEST_DIRNAME/common.bash
init_tests
# Setup parameters
SSFP_FILE="ssfp$EXT"
SSFP_SEQ="\"SSFP\": {
    \"TR\": 0.005,
    \"FA\": [12,12,30,30,45,45,65,65],
    \"PhaseInc\": [180,0,90,270,180,0,90,270]
}"
SIZE="16,16,16"
NOISE="0.001"
TOL="20"
cat > fm.json <<OUT
{
    $SSFP_SEQ
}
OUT
cat > simulate_fm.json <<OUT
{
    "PDFile": "PD$EXT",
    "T2File": "T2$EXT",
    "f0File": "f0$EXT",
    $SSFP_SEQ
}
OUT
}

@test "DESPOT2FM-Simulate" {
qinewimage PD$EXT --size "$SIZE" -f 1.0
[ -e PD$EXT ]
qinewimage T1$EXT --size "$SIZE" -f 1.0
[ -e T1$EXT ]
qinewimage T2$EXT --size "$SIZE" -g "0 0.025 0.125"
[ -e T2$EXT ]
qinewimage f0$EXT --size "$SIZE" -g "1 -150.0 150.0"
[ -e f0$EXT ]
qinewimage B1$EXT --size "$SIZE" -g "2 0.5 1.5"
[ -e B1$EXT ]
qidespot2fm --verbose --simulate=$NOISE --B1=B1$EXT T1$EXT $SSFP_FILE < simulate_fm.json
}

@test "DESPOT2FM-Fit" {
qidespot2fm --verbose --asym --B1=B1$EXT T1$EXT $SSFP_FILE < fm.json
qidiff --baseline=T2$EXT --input=FM_T2$EXT --noise=$NOISE --tolerance=$TOL --verbose
}
