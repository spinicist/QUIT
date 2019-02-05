# Copyright Tobias Wood 2017
# Tests for mcDESPOT

setup() {
load $BATS_TEST_DIRNAME/common.bash
init_tests
# Setup parameters
SPGR_FILE="spgr$EXT"
SSFP_FILE="ssfp$EXT"
SIZE="5,5,5"
SEQUENCE_GROUP="\
    \"Sequences\": [
        {
            \"SPGREcho\": {
                \"TR\": 0.0065,
                \"TE\": 0.003,
                \"FA\": [3,4,5,6,7,9,13,18]
            }
        },
        {
            \"SSFPEcho\": {
                \"TR\": 0.005,
                \"FA\": [12,16,21,27,33,40,51,68,12,16,21,27,33,40,51,68],
                \"PhaseInc\": [180,180,180,180,180,180,180,180,0,0,0,0,0,0,0,0]
            }
        }
    ]
"
NOISE="0.001"
TOL="300"
cat > 3pe.json <<OUT
{
    $SEQUENCE_GROUP
}
OUT
cat > simulate_3pe.json <<OUT
{
    "PDFile": "PD$EXT",
    "T1_mFile": "T1_m$EXT",
    "T2_mFile": "T2_m$EXT",
    "T1_ieFile": "T1_ie$EXT",
    "T2_ieFile": "T2_ie$EXT",
    "T1_csfFile": "T1_csf$EXT",
    "T2_csfFile": "T2_csf$EXT",
    "tau_mFile": "tau_m$EXT",
    "f_mFile": "f_m$EXT",
    "f_csfFile": "f_csf$EXT",
$SEQUENCE_GROUP
}
OUT
}

@test "mcDESPOT Three Pool w/ Echo-Time Correction Simulate" {
qinewimage --size "$SIZE" -f "1.0" PD$EXT
[ -e PD$EXT ]
qinewimage --size "$SIZE" -f "0.465" T1_m$EXT
[ -e T1_m$EXT ]
qinewimage --size "$SIZE" -f "0.026" T2_m$EXT
[ -e T2_m$EXT ]
qinewimage --size "$SIZE" -f "1.070" T1_ie$EXT
[ -e T1_ie$EXT ]
qinewimage --size "$SIZE" -f "0.117" T2_ie$EXT
[ -e T2_ie$EXT ]
qinewimage --size "$SIZE" -f "3.5" T1_csf$EXT
[ -e T1_csf$EXT ]
qinewimage --size "$SIZE" -f "1.0" T2_csf$EXT
[ -e T2_csf$EXT ]
qinewimage --size "$SIZE" -f "0.18" tau_m$EXT
[ -e tau_m$EXT ]
qinewimage --size "$SIZE" -g "0 0.05 0.25" f_m$EXT
[ -e f_m$EXT ]
qinewimage --size "$SIZE" -g "1 0.1  0.5" f_csf$EXT
[ -e f_csf$EXT ]
qinewimage --size "$SIZE" -g "2 -25. 25." f0$EXT
[ -e f0$EXT ]
qinewimage --size "$SIZE" -f "1" B1$EXT
[ -e B1$EXT ]
qimcdespot --verbose --simulate=$NOISE --model=3 --f0=f0$EXT --B1=B1$EXT $SPGR_FILE $SSFP_FILE < simulate_3p.json
}

@test "mcDESPOT Three Pool w/ Echo-Time Correction" {
qimcdespot --verbose --scale --f0=f0$EXT --model=3 --B1=B1$EXT $SPGR_FILE $SSFP_FILE < 3p.json
qidiff --baseline=f_m$EXT --input=3C_f_m$EXT --noise=$NOISE --tolerance=$TOL --verbose
qidiff --baseline=f_csf$EXT --input=3C_f_csf$EXT --noise=$NOISE --tolerance=$TOL --verbose
}
