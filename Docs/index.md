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

Pre-compiled binaries are provided for Linux and Mac OS X in a .tar.gz archive from http://github.com/spinicist/QUIT/releases.

Download the correct archive for your platform, untar it and then ensure that the binaries can be found via your `PATH` environment variable. Either edit your shell profile (e.g. `.bashrc`) and add the QUIT directory your `PATH` variable there, or copy the files to somewhere that will be on your path, e.g. `/usr/local/bin` (this will likely require `sudo` permissions).

The Linux binaries are built with Ubuntu 14.04 with GCC 6. If you need to run on an older version of Linux with a previous version of `glibc` then you will need to
compile from source.

##Compile From Source

See the [developer documentation](Developer.md)

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

By default, QUIT is compiled with support for NIFTI and NRRD formats. The preferred file-format is NIFTI for compatibility with FSL and SPM. By default QUIT will output `.nii.gz` files. This can be controlled by the `QUIT_EXT` environment variable. Valid values for this are any file extension supported by ITK that QUIT has been compiled to support, e.g. `.nii` or `.nrrd`, or the FSL values `NIFTI`, `NIFTI_PAIR`, `NIFTI_GZ`, `NIFTI_PAIR_GZ`.

The [ITK](http://itk.org) library supports a much wider variety of file formats, but adding support for all of these almost triples the size of the compiled binaries. Hence by default they are excluded. You can add support for more file formats by compiling QUIT yourself, see the [developer documentation](Developer.md). Note that ITK cannot write every format it can read (e.g. it can read Bruker 2dseq datasets, but it cannot write them).

## Scripting

The QUIT tools are designed to be used inside shell scripts or SGE jobs. 3 scripts are provided along with the binaries. `qi_composer.sh` uses the `qi_coil_combine` program to correctly combine complex data from multiple element RF coils. `qi_example_mcd_b1.sh` and `qi_example_mcd_hifi.sh` give example mcDESPOT processing pipelines. The first assumes that you have a good B1 map and have already registered all the input volumes to a common reference. The second shows one of doing this using FSL. Registration is beyond the scope of QUIT, so you must have additional tools installed for this.