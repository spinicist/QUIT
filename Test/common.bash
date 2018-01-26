#DATADIR="$BATS_TEST_DIRNAME/basic"
if [ -z $QI_TEST_DIR ]; then
    QI_TEST_DIR="$BATS_TMPDIR"
fi
BINDIR="$BATS_TEST_DIRNAME/.."
EXT=".nii.gz"
export QUIT_EXT="$EXT"

init_tests() {
    if [ -e ./qidiff ]; then
       # Looks like we are in the build dir, run locally
       echo "Found a local copy of QUIT, adding it to path"
       export PATH="$PWD/:$PATH"
    fi
    mkdir -p $QI_TEST_DIR
    cd $QI_TEST_DIR
    echo "Running tests in $QI_TEST_DIR"
    echo "Removing any previous image files"
    if ls *.nii* > /dev/null 2>&1; then
        rm *.nii*
    fi
}