![Logo](logo.png)

# Welcome

Welcome to the QUantitative Imaging Tools, a collection of C++ programs for analysing quantitative MR images. QUIT is divided into several modules, each of which contains numerous programs for processing images. The modules are:

* [Relaxometry](Relaxometry.md)
* [Perfusion](Perfusion.md)
* [Magnetization Transfer](MT.md)
* [Steady-State Free-Precession](SSFP.md)
* [Susceptibility](Susceptibility.md)
* [Stats / GLM](Stats.md)
* [Utilities](Utilities.md)

## Installation

**Binary Installation**

Compiled binaries are provided for Mac and Linux and can be found on the github releases page. These are provided in both tar (`.tar.gz`) and self-extracting tar (`.sh`) formats.

If you download the self-extracting tar version, the following command should install QUIT correctly:

```bash
sudo ./QUIT-2.0.0-Darwin.sh --prefix=/usr/local --exclude-subdir
```

If you download the `.tar.gz` version, you will need to copy the files to the correct locations. The usual location is `/usr/local`, however it is possible to locate them somewhere else. Ensure that the libraries (in /lib) are somewhere that the system can find them. On Mac, this means the directory with the libraries must be contained in the `DLYD_LIBRARY_PATH` environment variable.

**Compile From Source**

Compiling QUIT requires a recent version of CMake and a C++11 compatible compiler (GCC 5 or higher is recommended). If you have cloned the git repository and are compiling from source, it is recommended to run the `easy_build.sh` script the first time. This will ensure the external libraries required have been initialised properly and that the main QUIT build is pointed at them.

If you wish to customise your build, then the standard `mkdir build && cd build && ccmake ../` approach will work, but you will have to tell `cmake` where to find the external libraries (Args, Cereal, Eigen, Ceres & ITK). These are included as git submodules, but Ceres and ITK will require compilation first.

## General Usage

QUIT programs take their input as a combination of command-line arguments and text files passed to `stdin`. If you run a QUIT program with no arguments then will see a message like this:

```bash
~: qidespot1
SPGR FILE was not specified. Use --help to see usage.
>
```

If you then run it with `--help`, you will see some usage instructions:

```bash
~: qidespot1 --help
  qidespot1 {OPTIONS} [SPGR FILE]

    Calculates T1 maps from SPGR data
    http://github.com/spinicist/QUIT

  OPTIONS:

      SPGR FILE                         Path to SPGR data
      -h, --help                        Show this help message
      -v, --verbose                     Print more information
      -T[THREADS], --threads=[THREADS]  Use N threads (default=4, 0=hardware
                                        limit)
      -o[OUTPREFIX], --out=[OUTPREFIX]  Add a prefix to output filenames
      -b[B1], --B1=[B1]                 B1 map (ratio) file
      -m[MASK], --mask=[MASK]           Only process voxels within the mask
      -s[SUBREGION],
      --subregion=[SUBREGION]           Process subregion starting at voxel
                                        I,J,K with size SI,SJ,SK
      -r, --resids                      Write out residuals for each data-point
      -a[ALGO], --algo=[ALGO]           Choose algorithm (l/w/n)
      -i[ITERS], --its=[ITERS]          Max iterations for WLLS/NLLS (default
                                        15)
      -p[CLAMP PD], --clampPD=[CLAMP
      PD]                               Clamp PD between 0 and value
      -t[CLAMP T1], --clampT2=[CLAMP
      T1]                               Clamp T1 between 0 and value
      "--" can be used to terminate flag options and force all following
      arguments to be treated as positional options
```

The first line shows that DESPOT1 expects a single input image, in this case SPGR/FLASH/FFE, specified on the command line. However, if you naively run:

```bash
~: qidespot1 some_spgr_data.nii.gz
```

Nothing will happen. This is because most imaging formats do not store parameters that are required for quantitative imaging, e.g. the repetition time and flip-angles in their headers. QUIT programs hence expect to read this information in a text file passed to `stdin`. If you create a small text file `spgr.txt` containing the following:

```json
{
    "SPGR": {
        "TR": 0.01,
        "FA": [3, 18]
    }
}
```

and run the following:

```bash
~: qidespot1 some_spgr_data.nii.gz < spgr.txt
```

then (provided your input data does contain two volumes corresponding to flip-angles 3 and 18 degrees) then DESPOT1 will run, and you should see two files created (`D1_T1.nii.gz` and `D1_PD.nii.gz`). If you want to see what the programs are doing while running, specify the `--verbose` or `-v` options.

## File Formats

Because QUIT uses the [ITK](http://itk.org) library, it can read a wide variety of file-formats, including DICOM & NRRD. The preferred file-format is NIFTI for compatibility with FSL and SPM. By default QUIT will output `.nii.gz` files. This can be controlled by the `QUIT_EXT` environment variable. Valid values for this are any file extension that ITK can write, e.g. `.nii` or `.nrrd`, or the FSL values `NIFTI`, `NIFTI_PAIR`, `NIFTI_GZ`, `NIFTI_PAIR_GZ`. Note that ITK cannot write every format it can read (e.g. it can read Bruker 2dseq datasets, but it cannot write them).

## Scripting

The QUIT tools are designed to be used inside shell scripts or SGE jobs. 3 scripts are provided along with the binaries. `qi_composer.sh` uses the `qi_coil_combine` program to correctly combine complex data from multiple element RF coils. `qi_example_mcd_b1.sh` and `qi_example_mcd_hifi.sh` give example mcDESPOT processing pipelines. The first assumes that you have a good B1 map and have already registered all the input volumes to a common reference. The second shows one of doing this using FSL. Registration is beyond the scope of QUIT, so you must have additional tools installed for this.