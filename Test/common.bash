#DATADIR="$BATS_TEST_DIRNAME/basic"
DATADIR="$BATS_TMPDIR"
BINDIR="$BATS_TEST_DIRNAME/.."

init_tests() {
    if [ -e ./qidiff ]; then
       # Looks like we are in the build dir, run locally
       echo "Found a local copy of QUIT, adding it to path"
       export PATH="$PWD/:$PATH"
    fi
    mkdir -p $DATADIR
    cd $DATADIR
    echo "Running tests in $DATADIR"
    echo "Removing any previous image files"
    if ls *.nii* > /dev/null 2>&1; then
        rm *.nii*
    fi
}