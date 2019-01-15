# Copyright Tobias Wood 2017
# Tests for mcDESPOT

setup() {
load $BATS_TEST_DIRNAME/common.bash
init_tests
# Setup parameters
MTSAT_FILE="mtsat$EXT"
T1OBS_FILE="t1_obs$EXT"
MTSAT_SEQUENCE=''
NOISE="0.002"
TOL=40
cat >mtsat.json <<END_SIG
{
    "PDFile": "PD$EXT",
    "T1_fFile": "T1_f$EXT",
    "T2_fFile": "T2_f$EXT",
    "T2_bFile": "T2_b$EXT",
    "k_bfFile": "k_bf$EXT",
    "f_bFile": "f_b$EXT",
    "f0File": "f0$EXT",
    "B1File": "B1$EXT",
    "MTSat": {
        "TR": 0.055,
        "FA": 5,
        "sat_f0": [56360, 47180, 12060, 1000, 1000, 2750, 2770, 2790, 2890, 1000, 1000],
        "sat_angle": [332, 628, 628, 332, 333, 628, 628, 628, 628, 628, 628],
        "pulse": { "name": "Gauss", "Trf": 0.015, "p1": 0.416, "p2": 0.295 }
    }
}
END_SIG
}

@test "qMT Lineshapes" {
qi_lineshape --lineshape=Gaussian --frq_start=500 --frq_space=500 --frq_count=150 > gauss.json
[ -e gauss.json ]
qi_lineshape --lineshape=Lorentzian --frq_start=500 --frq_space=500 --frq_count=150 > lorentz.json
[ -e lorentz.json ]
qi_lineshape --lineshape=SuperLorentzian --frq_start=500 --frq_space=500 --frq_count=150 > superlorentz.json
[ -e superlorentz.json ]
}

@test "qMT Ramani Simulate" {
SIZE="10,10,10"
qinewimage --size "$SIZE" -f "1.0" PD$EXT
[ -e PD$EXT ]
qinewimage --size "$SIZE" -f "1.0" T1_f$EXT
[ -e T1_f$EXT ]
qinewimage --size "$SIZE" -f "0.03" T2_f$EXT
[ -e T2_f$EXT ]
qinewimage --size "$SIZE" -f "0.00001" T2_b$EXT
[ -e T2_b$EXT ]
qinewimage --size "$SIZE" -g "1 0.5 1.5" k_bf$EXT
[ -e k_bf$EXT ]
qinewimage --size "$SIZE" -g "0 0.05 0.25" f_b$EXT
[ -e f_b$EXT ]
qinewimage --size "$SIZE" -f "0.0" f0$EXT
[ -e f0$EXT ]
qinewimage --size "$SIZE" -g "2 0.8 1.1" B1$EXT
[ -e B1$EXT ]

qi_qmt --verbose --simulate=$NOISE --lineshape=superlorentz.json --file=mtsat.json $MTSAT_FILE $T1OBS_FILE
[ -e $MTSAT_FILE ]
[ -e $T1OBS_FILE ]
}

@test "qMT Ramani" {
qi_qmt --verbose --lineshape=superlorentz.json --file=mtsat.json $MTSAT_FILE $T1OBS_FILE
qidiff --verbose --baseline=f_b$EXT --input=QMT_f_b$EXT --noise=$NOISE --tolerance=$TOL
qidiff --verbose --baseline=k_bf$EXT --input=QMT_k_bf$EXT --noise=$NOISE --tolerance=$TOL
}