# Copyright Tobias Wood 2016
# Simple test scripts for QUIT programs

setup() {
    # DATADIR="$BATS_TMPDIR"
    DATADIR="$BATS_TEST_DIRNAME/d1"
    BINDIR="$BATS_TEST_DIRNAME/.."
    # Setup parameters
    SPGR_FILE="spgr.nii"
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
    command mkdir -p $DATADIR
    command cd $DATADIR
    if [ "$(ls -A ./)" ]; then
        rm *.nii *.txt
    fi

    SIZE="32 32 32"
    command $BINDIR/qinewimage --size "$SIZE" -g "1 0.8 1.0" PD.nii
    [ -e PD.nii ]
    command $BINDIR/qinewimage --size "$SIZE" -g "0 0.5 1.5" T1.nii
    [ -e T1.nii ]
    command $BINDIR/qinewimage --size "$SIZE" -f "0.05" T2.nii
    [ -e T2.nii ]
    command $BINDIR/qinewimage --size "$SIZE" -f "0.0" f0.nii
    [ -e f0.nii ]
    command $BINDIR/qinewimage --size "$SIZE" -g "2 0.75 1.25" B1.nii
    [ -e B1.nii ]
}

@test "Generate test signals" {
    command cd $DATADIR
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
    [ -e $SPGR_FILE ]
    [ -e $MPRAGE_FILE ]
    [ -e $AFI_FILE ]
}

@test "Process HIFI" {
    command cd $DATADIR
    echo \
"$SPGR_FLIP
$SPGR_TR
$MPRAGE_FLIP
$SPGR_TR
$MPRAGE_SEGSIZE
$MPRAGE_KZERO
$MPRAGE_INV
$MPRAGE_DELAY
$MPRAGE_EFF" > hifi.txt
    command $BINDIR/qidespot1hifi $SPGR_FILE $MPRAGE_FILE -M -n -T1 -v < hifi.txt
}

@test "Process AFI" {
    command cd $DATADIR
    echo \
"$SPGR_FLIP
$SPGR_TR" > despot1.txt
    command $BINDIR/qiafi $AFI_FILE -v
    command $BINDIR/qidespot1 $SPGR_FILE -v -n -bAFI_B1.nii -o AFI -an < despot1.txt
}
