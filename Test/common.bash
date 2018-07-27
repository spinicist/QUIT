echo "$PATH"

if [ -z $QI_TEST_DIR ]; then
    QI_TEST_DIR="$BATS_TMPDIR"
fi
EXT=".nii.gz"
export QUIT_EXT="$EXT"

init_tests() {
    mkdir -p $QI_TEST_DIR
    cd $QI_TEST_DIR
    echo "Running tests in $QI_TEST_DIR"
    # echo "Removing any previous image files"
    # if ls *.nii* > /dev/null 2>&1; then
    #     rm *.nii*
    # fi
}