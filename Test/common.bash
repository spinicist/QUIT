#DATADIR="$BATS_TEST_DIRNAME/basic"
if [ -z $QI_TEST_DIR ]; then
    QI_TEST_DIR="$BATS_TMPDIR"
fi
EXT=".nii.gz"
export QUIT_EXT="$EXT"

if [[ $OSTYPE == darwin* ]]; then
    # Work around SIP removing DYLD_LIBRARY_PATH in child processes
    if [ -e $HOME/.bashrc ]; then
        source $HOME/.bashrc
    elif [ -e $HOME/.bash_profile ]; then
        source $HOME/.bash_profile
    fi
fi

init_tests() {
    mkdir -p $QI_TEST_DIR
    cd $QI_TEST_DIR
    echo "Running tests in $QI_TEST_DIR"
    echo "Removing any previous image files"
    if ls *.nii* > /dev/null 2>&1; then
        rm *.nii*
    fi
}