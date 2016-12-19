# Copyright Tobias Wood 2016
# Simple test scripts for QUIT programs

setup() {
DATADIR="d1"
BINDIR=$PWD/..
mkdir -p $DATADIR
cd $DATADIR
# Setup parameters
SPGR_FILE="spgr2.nii"
SPGR_FLIP="3 3 20 20"
SPGR_TR="0.01"
MPRAGE_FILE="mprage.nii"
MPRAGE_FLIP="5"
MPRAGE_SEGSIZE="64"
MPRAGE_KZERO="0"
MPRAGE_INV="0.45"
MPRAGE_DELAY="0.0"
MPRAGE_EFF="0.95"
AFI_FILE="afi.nii"
AFI_PAR="55.
0.02 0.1"
}

@test "Generate test parameter images" {
    if [ "$(ls -A ./)" ]; then
        rm *
    fi

    SIZE="32 32 32"
    run $BINDIR/qinewimage --size "$SIZE" -g "1 0.8 1.0" PD.nii
    [ "$status" -eq 0 ]
    [ -e PD.nii ]
    run $BINDIR/qinewimage --size "$SIZE" -g "0 0.5 1.5" T1.nii
    [ "$status" -eq 0 ]
    [ -e T1.nii ]
    run $BINDIR/qinewimage --size "$SIZE" -f "0.05" T2.nii
    [ "$status" -eq 0 ]
    [ -e T2.nii ]
    run $BINDIR/qinewimage --size "$SIZE" -f "0.0" f0.nii
    [ "$status" -eq 0 ]
    [ -e f0.nii ]
    run $BINDIR/qinewimage --size "$SIZE" -g "2 0.75 1.25" B1.nii
    [ "$status" -eq 0 ]
    [ -e B1.nii ]
}

@test "Generate test signals" {
echo \
"PD.nii
T1.nii
T2.nii
f0.nii
B1.nii
$SPGR_FILE
SPGR
$SPGR_FLIP
$SPGR_TR
$MPRAGE_FILE
MPRAGE
$MPRAGE_FLIP
$SPGR_TR
$MPRAGE_SEGSIZE
$MPRAGE_KZERO
$MPRAGE_INV
$MPRAGE_DELAY
$MPRAGE_EFF
$AFI_FILE
AFI
$AFI_PAR
END" > signal_input.txt

    NOISE="0.002"
    command $BINDIR/qisignal --model=1 -v -n --noise=$NOISE < signal_input.txt
    [ -e $SPGR2_FILE ]
    [ -e $SPGR4_FILE ]
    [ -e $MPRAGE_FILE ]
    [ -e $AFI_FILE ]
}
