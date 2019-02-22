# Copyright Tobias Wood 2017
# Tests for perfusion (ASL)

setup() {
load $BATS_TEST_DIRNAME/common.bash
init_tests
}

@test "Perfusion (ASL)" {

SIZE="3,3,3,2"
qinewimage --verbose --dims=4 --size="${SIZE},2" --step="3 1 1.06 2" asl$EXT
qinewimage --verbose --dims=3 --size="${SIZE}" --fill=147.355 ref_cbf$EXT
qi_asl --verbose asl$EXT <<END_INPUT
{
    "CASL" : {
        "TR": 4.0,
        "label_time": 3.0,
        "post_label_delay": [0.3]
    }
}
END_INPUT
qidiff --baseline=ref_cbf$EXT --input=CASL_CBF$EXT --noise=0 --tolerance=1.e-3 --verbose
}

@test "Perfusion (ASE)" {
SIZE="9,9,9"
qinewimage S0$EXT --size="$SIZE" -f 100
qinewimage dT$EXT --size="$SIZE" -g "0 -0.025 0.025"
qinewimage R2p$EXT --size="$SIZE" -g "1 1.0 3.0"
qinewimage DBV$EXT --size="$SIZE" -g "2 0.005 0.025"

cat > input.json <<END
{
    "MultiEchoFlex" : {
        "TR" : 2.5,
        "TE" : [ -0.05, 
                 -0.045,
                 -0.04,
                 -0.035,
                 -0.03,
                 -0.025,
                 -0.02,
                 -0.015,
                 -0.01,
                 -0.005,
                  0.0,
                  0.005,
                  0.01,
                  0.015,
                  0.02,
                  0.025,
                  0.03,
                  0.035,
                  0.04,
                  0.045,
                  0.05 ]
    },
    "S0File" : "S0$EXT",
    "dTFile" : "dT$EXT",
    "R2pFile" : "R2p$EXT",
    "DBVFile" : "DBV$EXT"
}
END
SPIN_FILE="me$EXT"
NOISE="0.01"
qi_ase_oef --verbose --simulate=$NOISE $SPIN_FILE --threads=1 < input.json
qi_ase_oef --verbose $SPIN_FILE --threads=1 < input.json
qidiff --baseline=R2p$EXT --input=ASE_R2p$EXT --noise=$NOISE --tolerance=5 --verbose
qidiff --baseline=DBV$EXT --input=ASE_DBV$EXT --noise=$NOISE --tolerance=10 --verbose
}

@test "Perfusion (Z-Shim)" {
# Pythagoras 3,4,5 triangle!
qinewimage zshim_in$EXT --dims=4 --size="2,2,2,2" --step="3 3 4 2"
qinewimage zshim_ref$EXT --size="2,2,2" --fill="5"
qi_zshim --verbose --zshims=2 zshim_in$EXT
qidiff --baseline=zshim_ref$EXT --input=zshim_in_zshim$EXT --abs --tolerance=0.1 --verbose
}