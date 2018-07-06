# Copyright Tobias Wood 2017
# Simple test scripts for QUIT programs

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "DESPOT1-HIFI" {

# Setup parameters
SPGR_FILE="spgr$EXT"
MPRAGE_FILE="mprage$EXT"
SPGR=' "SPGR": { "TR": 0.01, "FA": [3, 20] }'
echo $SPGR
MPRAGE=' "MPRAGE": { "FA": 5, "TR": 0.01, "TI": 0.45, "TD": 0, "eta": 1, "ETL": 64, "k0": 0 }'

SIZE="16,16,16"
NOISE="0.01"
qinewimage --size "$SIZE" -g "1 0.8 1.0" PD$EXT
qinewimage --size "$SIZE" -g "0 0.5 1.5" T1$EXT
qisignal --model=1 -v --noise=$NOISE $SPGR_FILE $MPRAGE_FILE << SIGNAL
{
    "PD": "PD$EXT",
    "T1": "T1$EXT",
    "T2": "",
    "f0": "",
    "B1": "",
    "Sequences": [
        { $SPGR },
        { $MPRAGE }
    ]
}
SIGNAL

qidespot1hifi --verbose $SPGR_FILE $MPRAGE_FILE <<HIFI
{
    $SPGR, $MPRAGE
}
HIFI

qidiff --baseline=T1.nii --input=HIFI_T1.nii --noise=$NOISE --tolerance=40 --verbose

}
