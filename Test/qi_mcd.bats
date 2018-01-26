# Copyright Tobias Wood 2017
# Tests for mcDESPOT

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "2C mcDESPOT" {

SIZE="4,4,4"
qinewimage --size "$SIZE" -f "1.0" PD$EXT
qinewimage --size "$SIZE" -f "0.465" T1_m$EXT
qinewimage --size "$SIZE" -f "0.026" T2_m$EXT
qinewimage --size "$SIZE" -f "1.070" T1_ie$EXT
qinewimage --size "$SIZE" -f "0.117" T2_ie$EXT
qinewimage --size "$SIZE" -f "0.18" tau_m$EXT
qinewimage --size "$SIZE" -g "0 0. 200." f0$EXT
qinewimage --size "$SIZE" -g "1 0.75 1.25" B1$EXT
qinewimage --size "$SIZE" -g "2 0.05 0.25" f_m$EXT

# Setup parameters
SPGR_FILE="spgr$EXT"
SPGR_FLIP="3,4,5,6,7,9,13,18"
SPGR_TR="0.0065"
SSFP_FILE="ssfp$EXT"
SSFP_FLIP="12,16,21,27,33,40,51,68,12,16,21,27,33,40,51,68"
SSFP_PINC="180,180,180,180,180,180,180,180,0,0,0,0,0,0,0,0"
SSFP_TR="0.005"

SEQUENCE_GROUP="\
    \"SequenceGroup\": {
        \"sequences\": [
            {
                \"SPGR\": {
                    \"TR\": $SPGR_TR,
                    \"FA\": [$SPGR_FLIP]
                }
            },
            {
                \"SSFP\": {
                    \"TR\": $SSFP_TR,
                    \"FA\": [$SSFP_FLIP],
                    \"PhaseInc\": [$SSFP_PINC]
                }
            }
        ]
    }
"

NOISE="0.002"
qisignal --model=2 -v --noise=$NOISE $SPGR_FILE $SSFP_FILE << END_SIG
{
    "PD": "PD$EXT",
    "T1_m": "T1_m$EXT",
    "T2_m": "T2_m$EXT",
    "T1_ie": "T1_ie$EXT",
    "T2_ie": "T2_ie$EXT",
    "tau_m": "tau_m$EXT",
    "f_m": "f_m$EXT",
    "f0": "f0$EXT",
    "B1": "B1$EXT",
$SEQUENCE_GROUP
}
END_SIG

qimcdespot $OPTS -M2 -bB1$EXT -ff0$EXT -r -v $SPGR_FILE $SSFP_FILE << END_MCD
{
$SEQUENCE_GROUP
}
END_MCD
qidiff --baseline=f_m$EXT --input=2C_f_m$EXT --noise=$NOISE --tolerance=250 --verbose

}

@test "3C mcDESPOT" {

SIZE="4,4,4"
qinewimage --size "$SIZE" -f "1.0" PD$EXT
qinewimage --size "$SIZE" -f "0.465" T1_m$EXT
qinewimage --size "$SIZE" -f "0.026" T2_m$EXT
qinewimage --size "$SIZE" -f "1.070" T1_ie$EXT
qinewimage --size "$SIZE" -f "0.117" T2_ie$EXT
qinewimage --size "$SIZE" -f 4.0 T1_csf$EXT
qinewimage --size "$SIZE" -f 2.5 T2_csf$EXT
qinewimage --size "$SIZE" -f "0.18" tau_m$EXT
qinewimage --size "$SIZE" -g "0 0. 200." f0$EXT
qinewimage --size "$SIZE" -g "1 0.75 1.25" B1$EXT
qinewimage --size "$SIZE" -g "2 0.05 0.25" f_m$EXT
qinewimage --size "$SIZE" -f 0.05 f_csf$EXT

# Setup parameters
SPGR_FILE="spgr$EXT"
SPGR_FLIP="3,4,5,6,7,9,13,18"
SPGR_TR="0.0065"
SSFP_FILE="ssfp$EXT"
SSFP_FLIP="12,16,21,27,33,40,51,68,12,16,21,27,33,40,51,68"
SSFP_PINC="180,180,180,180,180,180,180,180,0,0,0,0,0,0,0,0"
SSFP_TR="0.005"

SEQUENCE_GROUP="\
    \"SequenceGroup\": {
        \"sequences\": [
            {
                \"SPGR\": {
                    \"TR\": $SPGR_TR,
                    \"FA\": [$SPGR_FLIP]
                }
            },
            {
                \"SSFP\": {
                    \"TR\": $SSFP_TR,
                    \"FA\": [$SSFP_FLIP],
                    \"PhaseInc\": [$SSFP_PINC]
                }
            }
        ]
    }
"

echo "$INPUT" > mcd_input.txt
NOISE="0.002"
qisignal --model=3 -v --noise=$NOISE $SPGR_FILE $SSFP_FILE << END_SIG
{
    "PD" : "PD$EXT",
    "T1_m" : "T1_m$EXT",
    "T2_m" : "T2_m$EXT",
    "T1_ie" : "T1_ie$EXT",
    "T2_ie" : "T2_ie$EXT",
    "T1_csf" : "T1_csf$EXT",
    "T2_csf" : "T2_csf$EXT",
    "tau_m" : "tau_m$EXT",
    "f_m" : "f_m$EXT",
    "f_csf" : "f_csf$EXT",
    "f0" : "f0$EXT",
    "B1" : "B1$EXT",
$SEQUENCE_GROUP
}
END_SIG

qimcdespot $OPTS -M3 -bB1$EXT -ff0$EXT -r -v $SPGR_FILE $SSFP_FILE << END_MCD
{ $SEQUENCE_GROUP }
END_MCD
qidiff --baseline=f_m$EXT --input=3C_f_m$EXT --noise=$NOISE --tolerance=250 --verbose

}