# Copyright Tobias Wood 2017
# Tests for perfusion (ASL)

setup() {
load $BATS_TEST_DIRNAME/common.bash
init_tests
}

@test "Perfusion (ASL)" {

SIZE="32,32,32,2"
qinewimage --verbose --dims=4 --size="$SIZE" --step="3 1 1.06 2" asl$EXT
qi_asl --verbose asl$EXT <<END_INPUT
{
    "CASL" : {
        "TR": 4.0,
        "label_time": 3.0,
        "post_label_delay": [0.3]
    }
}
END_INPUT
}

@test "Perfusion (ASE)" {
# Currently, this test uses a T2 decay curve as the input and checks that R2' is equal to 1/T2
# and DBV is equal to zero
SIZE="16,16,16"
qinewimage PD$EXT --size="$SIZE" -f 1
qinewimage T1$EXT --size="$SIZE" -f 1
qinewimage T2$EXT --size="$SIZE" -f 0.1
qinewimage R2$EXT --size="$SIZE" -f 10
qinewimage DBV$EXT --size="$SIZE" -f 0.0

SPIN_SEQUENCE="
\"MultiEchoFlex\" : {
    \"TR\" : 2.5,
    \"TE\" : [ 0.00, 0.00, 0.00, 0.05, 0.1, 0.15, 0.3, 0.45, 0.6 ]
}
"
SPIN_FILE="me$EXT"
NOISE="0.01"
qimultiecho --verbose --simulate=$NOISE $SPIN_FILE << END_SIG
{
    "PDFile" : "PD$EXT",
    "T2File" : "T2$EXT",
    $SPIN_SEQUENCE
}
END_SIG

qi_ase_oef --verbose $SPIN_FILE --threads=1 << END_INPUT
{
    $SPIN_SEQUENCE
}
END_INPUT
qidiff --baseline=R2$EXT --input=me_R2prime$EXT --noise=$NOISE --tolerance=5 --verbose
qidiff --baseline=DBV$EXT --input=me_DBV$EXT --noise=$NOISE --tolerance=250 --abs --verbose
}

@test "Perfusion (Z-Shim)" {
# Pythagoras 3,4,5 triangle!
qinewimage zshim_in$EXT --dims=4 --size="2,2,2,2" --step="3 3 4 2"
qinewimage zshim_ref$EXT --size="2,2,2" --fill="5"
qi_zshim --verbose --zshims=2 zshim_in$EXT
qidiff --baseline=zshim_ref$EXT --input=zshim_in_zshim$EXT --abs --tolerance=0.1 --verbose
}