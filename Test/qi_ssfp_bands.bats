# Copyright Tobias Wood 2018
# Simple test scripts for QUIT programs

setup() {
load $BATS_TEST_DIRNAME/common.bash
init_tests
# Setup parameters
SIZE="8,8,8"
NOISE="0.005"
TOL="10"
SSFP_SEQUENCE="
\"SSFP\": {
            \"TR\": 0.01,
            \"FA\": [30,30,30,30],
            \"PhaseInc\": [180,0,90,270]
}"
}

@test "SSFP-Ellipse Simulate" {
# Borrow qi_ssfp_ellipse
qinewimage --size "$SIZE" -f 1.0 G$EXT
[ -e G$EXT ]
qinewimage --size "$SIZE" -g "0 0.9 0.99" a$EXT
[ -e a$EXT ]
qinewimage --size "$SIZE" -g "1 0.3 0.8" b$EXT
[ -e b$EXT ]
qinewimage --size "$SIZE" -g "2 -3 3" theta_0$EXT
[ -e theta_0$EXT ]
qinewimage --size "$SIZE" -f 0.0 phi_rf$EXT
[ -e phi_rf$EXT ]
cat > simulate_ellipse.json <<OUT
{
    "GFile": "G$EXT",
    "aFile": "a$EXT",
    "bFile": "b$EXT",
    "theta_0File": "theta_0$EXT",
    "phi_rfFile": "phi_rf$EXT",
    $SSFP_SEQUENCE
}
OUT
qi_ssfp_ellipse --verbose --simulate=$NOISE ssfp$EXT < simulate_ellipse.json
[ -e ssfp$EXT ]
qicomplex -x ssfp$EXT -M ssfp_mag$EXT
[ -e ssfp_mag$EXT ]
}

@test "SSFP Mag Banding" {
MAG_TOL="300"
qi_ssfp_bands ssfp_mag$EXT --method=M --magnitude --verbose
qidiff --baseline=G$EXT --input=ssfp_mag_Max$EXT --noise=$NOISE --tolerance=$MAG_TOL --verbose
qi_ssfp_bands ssfp_mag$EXT --method=N --magnitude --verbose
qidiff --baseline=G$EXT --input=ssfp_mag_MagMean$EXT --noise=$NOISE --tolerance=$MAG_TOL --verbose
}

@test "SSFP GS Banding" {
    qi_ssfp_bands ssfp$EXT --method=G --magnitude --alt-order --verbose
    qidiff --baseline=G$EXT --input=ssfp_GSL$EXT --noise=$NOISE --tolerance=$TOL --verbose
    qi_ssfp_bands ssfp$EXT --method=G -2 --magnitude --alt-order --verbose
    qidiff --baseline=G$EXT --input=ssfp_GSL2$EXT --noise=$NOISE --tolerance=$TOL --verbose
}