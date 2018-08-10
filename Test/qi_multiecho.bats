# Copyright Tobias Wood 2017
# Tests for relaxometry programs

setup() {
load $BATS_TEST_DIRNAME/common.bash
init_tests
# Setup parameters
SIZE="16,16,16"
SPIN_FILE="me$EXT"
SPIN_SEQ="\"MultiEcho\" : { \"TR\": 2.5, \"TE1\": 0.005, \"ESP\": 0.005, \"ETL\": 16 }"
# Abuse the MultiEchoFlex sequence to get a two-volume output
SPIN_SEQ2="\"MultiEchoFlex\" : { \"TR\": 2.5, \"TE\": [0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.075,0.08,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.075,0.08] }"
NOISE="0.001"
TOL="25"
cat > me.json <<OUT
{
    $SPIN_SEQ
}
OUT
cat > simulate_me.json <<OUT
{
    "PDFile": "PD$EXT",
    "T2File": "T2$EXT",
    $SPIN_SEQ
}
OUT
cat > simulate_me2.json <<OUT
{
    "PDFile": "PD$EXT",
    "T2File": "T2$EXT",
    $SPIN_SEQ2
}
OUT
}

@test "MultiEcho-Simulate" {
qinewimage PD$EXT --size "$SIZE" -f 1
[ -e PD$EXT ]
qinewimage T2$EXT --size "$SIZE" -g "1 0.01 0.1"
[ -e T2$EXT ]
qimultiecho --verbose --simulate=$NOISE $SPIN_FILE < simulate_me.json
[ -e $SPIN_FILE ]
}

@test "MultiEcho-LLS" {
qimultiecho --verbose --resids $SPIN_FILE < me.json
qidiff --baseline=T2$EXT --input=ME_T2$EXT --noise=$NOISE --tolerance=$TOL --verbose
}

@test "MultiEcho-ARLO" {
qimultiecho --verbose --resids --algo=a $SPIN_FILE < me.json
qidiff --baseline=T2$EXT --input=ME_T2$EXT --noise=$NOISE --tolerance=$TOL --verbose
}

@test "MultiEcho-NLLS" {
qimultiecho --verbose --resids --algo=n $SPIN_FILE < me.json
qidiff --baseline=T2$EXT --input=ME_T2$EXT --noise=$NOISE --tolerance=$TOL --verbose
}

@test "MultiEcho-MultiVolume-Simulate" {
qinewimage PD$EXT --size "$SIZE" -f 1
[ -e PD$EXT ]
qinewimage T2$EXT --size "$SIZE" -g "1 0.01 0.1"
[ -e T2$EXT ]
qimultiecho --verbose --simulate=$NOISE me2$EXT < simulate_me2.json
[ -e $SPIN_FILE ]
}

@test "MultiEcho-MultiVolume-LLS" {
qimultiecho --verbose --resids me2$EXT < me.json
qidiff --baseline=T2$EXT --input=ME_T2$EXT --noise=$NOISE --tolerance=$TOL --verbose
[ $( qihdr --size=4 ME_T2$EXT) -eq "2" ]
}

@test "MultiEcho-MultiVolume-ARLO" {
qimultiecho --verbose --resids --algo=a me2$EXT < me.json
qidiff --baseline=T2$EXT --input=ME_T2$EXT --noise=$NOISE --tolerance=$TOL --verbose
[ $( qihdr --size=4 ME_T2$EXT) -eq "2" ]
}

@test "MultiEcho-MultiVolume-NLLS" {
qimultiecho --verbose --resids --algo=n me2$EXT < me.json
qidiff --baseline=T2$EXT --input=ME_T2$EXT --noise=$NOISE --tolerance=$TOL --verbose
[ $( qihdr --size=4 ME_T2$EXT) -eq "2" ]
}