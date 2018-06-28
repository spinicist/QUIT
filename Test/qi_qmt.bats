# Copyright Tobias Wood 2017
# Tests for mcDESPOT

setup() {
    load $BATS_TEST_DIRNAME/common.bash
    init_tests
}

@test "qMT Ramani" {
SIZE="3,3,3"
qinewimage --size "$SIZE" -f "1.0" PD$EXT
qinewimage --size "$SIZE" -f "0.1" T1_f$EXT
qinewimage --size "$SIZE" -f "0.00001" T2_f$EXT
qinewimage --size "$SIZE" -f "1.00" T1_b$EXT
qinewimage --size "$SIZE" -f "0.08" T2_b$EXT
qinewimage --size "$SIZE" -f "6" k_bf$EXT
qinewimage --size "$SIZE" -g "0 0.05 0.25" f_b$EXT
qinewimage --size "$SIZE" -g "1 0. 200." f0$EXT
qinewimage --size "$SIZE" -g "2 0.75 1.25" B1$EXT


# Setup parameters
MTSAT_FILE="mtsat$EXT"
MTSAT_SEQUENCE='{ "MTSat": {
    "TR": 30,
    "FA": 5,
    "sat_f0": [1000, 5000, 10000, 1000, 5000, 10000],
    "sat_angle": [360, 360, 360, 720, 720, 720],
    "pulse": { "name": "Gauss", "FAnom": 90, "Trf": 15, "intB1": 0.005904, "intB1sq": 0.005904 }
}}'
NOISE="0.002"

cat >mtsat.json <<END_SIG
{
    "PD": "PD$EXT",
    "T1_f": "T1_f$EXT",
    "T2_f": "T2_f$EXT",
    "T1_b": "T1_b$EXT",
    "T2_b": "T2_b$EXT",
    "k_bf": "k_bf$EXT",
    "f_b": "f_b$EXT",
    "f0": "f0$EXT",
    "B1": "B1$EXT",
    "SequenceGroup": {
        "sequences": [ $MTSAT_SEQUENCE ]
    }
}
END_SIG

qisignal --model=Ramani -v -T1 --noise=$NOISE $MTSAT_FILE << END_SIG
{
    "PD": "PD$EXT",
    "T1_f": "T1_f$EXT",
    "T2_f": "T2_f$EXT",
    "T1_b": "T1_b$EXT",
    "T2_b": "T2_b$EXT",
    "k_bf": "k_bf$EXT",
    "f_b": "f_b$EXT",
    "f0": "f0$EXT",
    "B1": "B1$EXT",
    "lineshape" : { "Gaussian" : {} },
    "SequenceGroup": {
        "sequences": [ $MTSAT_SEQUENCE ]
    }
}
END_SIG
return 1
}