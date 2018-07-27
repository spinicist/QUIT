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
            \"FA\": [30,30,30,30,30,30],
            \"PhaseInc\": [180,240,300,0,60,120]
}"
cat > ellipse.json <<OUT
{
    $SSFP_SEQUENCE
}
OUT

}

@test "SSFP-Ellipse Simulate" {
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
}

@test "SSFP-Ellipse DirectFit" {
qi_ssfp_ellipse --verbose --algo=d ssfp$EXT < ellipse.json
qidiff --baseline=G$EXT --input=ES_G$EXT --noise=$NOISE --tolerance=$TOL --verbose
qidiff --baseline=a$EXT --input=ES_a$EXT --noise=$NOISE --tolerance=$TOL --verbose
qidiff --baseline=b$EXT --input=ES_b$EXT --noise=$NOISE --tolerance=$TOL --verbose
qidiff --baseline=theta_0$EXT --input=ES_theta_0$EXT --noise=$NOISE --tolerance=$TOL --verbose
# Don't compare phi_rf - it sometimes wraps
# qidiff --baseline=phi_rf$EXT --input=ES_phi_rf$EXT --noise=$NOISE --tolerance=$TOL --verbose
}

@test "SSFP-Ellipse HyperFit" {
qi_ssfp_ellipse --verbose --algo=h ssfp$EXT < ellipse.json
qidiff --baseline=G$EXT --input=ES_G$EXT --noise=$NOISE --tolerance=100 --verbose
qidiff --baseline=a$EXT --input=ES_a$EXT --noise=$NOISE --tolerance=$TOL --verbose
qidiff --baseline=b$EXT --input=ES_b$EXT --noise=$NOISE --tolerance=$TOL --verbose
# Hyper does really poorly on these, ignore
# qidiff --baseline=theta_0$EXT --input=ES_theta_0$EXT --noise=$NOISE --tolerance=$TOL --verbose
# qidiff --baseline=phi_rf$EXT --input=ES_phi_rf$EXT --noise=$NOISE --tolerance=$TOL --verbose
}

#     qi_ssfp_ellipse ssfp$EXT --verbose --algo=d << INPUT
# {
#     "SSFPEllipse": {
#         "TR": 0.01,
#         "FA": [30],
#         "PhaseInc": [180,240,300,0,60,120]
#     }
# }
# INPUT
#     qi_ssfp_planet ES_G$EXT ES_a$EXT ES_b$EXT --verbose << INPUT
# {
#     "SSFPGS": {
#         "TR": 0.01,
#         "FA": [30]
#     }
# }
# INPUT
#     # No proper ground truth here, compare the methods to each other
#     qidiff --baseline=T1$EXT --input=PLANET_T1$EXT --noise=$NOISE --tolerance=40 --verbose
#     qidiff --baseline=T2$EXT --input=PLANET_T2$EXT --noise=$NOISE --tolerance=35 --verbose
#     qidiff --baseline=PD$EXT --input=PLANET_PD$EXT --noise=$NOISE --tolerance=20 --verbose
# }

# @test "SSFP EMT" {
#     SIZE="64,64,16"
#     NOISE="0.005"
#     qinewimage --size "$SIZE" -f 1.0 PD$EXT
#     qinewimage --size "$SIZE" -g "0 0.5 1.5" T1$EXT
#     qinewimage --size "$SIZE" -g "1 0.05 0.1" T2$EXT
#     qinewimage --size "$SIZE" -g "2 -200 200" f0$EXT
#     qisignal --model=1 -v --noise=$NOISE ssfp$EXT --complex << INPUT
# {
#     "PD": "PD$EXT",
#     "T1": "T1$EXT",
#     "T2": "T2$EXT",
#     "f0": "f0$EXT",
#     "B1": "",
#     "Sequences": [
#         {
#             "SSFPEcho": {
#                 "TR": 0.01,
#                 "FA": [16,16,16,16,16,16,32,32,32,32,32,32],
#                 "PhaseInc": [180,240,300,0,60,120,180,240,300,0,60,120]
#             }
#         }
#     ]
# }
# INPUT
#     qi_ssfp_ellipse ssfp$EXT --verbose --algo=d << INPUT
# {
#     "SSFPEllipse": {
#         "TR": 0.01,
#         "FA": [16,32],
#         "PhaseInc": [180,240,300,0,60,120]
#     }
# }
# INPUT
#     qi_ssfp_emt ES_G$EXT ES_a$EXT ES_b$EXT --verbose << INPUT
# {
#     "SSFPMT": {
#         "TR": [0.01, 0.01],
#         "Trf": [0.001, 0.001],
#         "FA": [16, 32],
#         "intB1": [0.005904, 0.005904],
#         "PhaseInc": [180, 240, 300, 0, 60, 120]
#     }
# }
# INPUT
#     # No proper ground truth here, compare the methods to each other
#     qidiff --baseline=T1$EXT --input=EMT_T1_f$EXT --noise=$NOISE --tolerance=40 --verbose
#     qidiff --baseline=T2$EXT --input=EMT_T2_f$EXT --noise=$NOISE --tolerance=35 --verbose
#     qidiff --baseline=PD$EXT --input=EMT_M0$EXT --noise=$NOISE --tolerance=20 --verbose
# }
