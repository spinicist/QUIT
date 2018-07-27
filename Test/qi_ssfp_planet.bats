# Copyright Tobias Wood 2018
# Simple test scripts for QUIT programs

setup() {
load $BATS_TEST_DIRNAME/common.bash
init_tests
# Setup parameters
SIZE="8,8,8"
NOISE="0.005"
TOL="30"
SSFP_SEQUENCE="
\"SSFP\": {
            \"TR\": 0.01,
            \"FA\": [10, 20, 30],
            \"PhaseInc\": [180, 180, 180]
}"
cat > planet.json <<OUT
{
    $SSFP_SEQUENCE
}
OUT

}

@test "SSFP-PLANET Simulate" {
qinewimage --size "$SIZE" -g "0 0.9 1.1" PD$EXT
[ -e PD$EXT ]
qinewimage --size "$SIZE" -g "1 0.9 1.5" T1$EXT
[ -e T1$EXT ]
qinewimage --size "$SIZE" -g "2 0.02 0.1" T2$EXT
[ -e T2$EXT ]
cat > simulate_planet.json <<OUT
{
    "PDFile": "PD$EXT",
    "T1File": "T1$EXT",
    "T2File": "T2$EXT",
    $SSFP_SEQUENCE
}
OUT
qi_ssfp_planet --verbose --simulate=$NOISE G$EXT a$EXT b$EXT < simulate_planet.json
[ -e G$EXT ]
[ -e a$EXT ]
[ -e b$EXT ]
}

@test "SSFP-PLANET Fit" {
qi_ssfp_planet --verbose G$EXT a$EXT b$EXT < planet.json
qidiff --baseline=PD$EXT --input=PLANET_PD$EXT --noise=$NOISE --tolerance=$TOL --verbose
qidiff --baseline=T1$EXT --input=PLANET_T1$EXT --noise=$NOISE --tolerance=$TOL --verbose
qidiff --baseline=T2$EXT --input=PLANET_T2$EXT --noise=$NOISE --tolerance=$TOL --verbose
}
