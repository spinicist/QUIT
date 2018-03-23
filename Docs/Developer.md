# Developer

Below are a few introductory notes for anyone who wants to compile QUIT from source or might be interested in contributing to QUIT. Although I hope that most QUIT code is fairly clear, the following is provided to give a high-level overview of how most QUIT programs are structured.

## Basic Requirements

1. A C++11 compliant compiler (GCC 4.8.0+, Clang 3.1+)
2. CMake version 3.2 or higher (http://www.cmake.org)

WARNING - You will require a recent compiler for DESPOT as it uses C++11 features. GCC 4.8.0 is a MINIMUM as earlier versions do not support C++11 threads - however you will not be able to build the `Stats` module without GCC 5.0.0 or later. Clang 3.1 will also work. You can find out the version of your system gcc by running `gcc -v`.

## External Libraries

QUIT is built using several C++ libraries. These are currently included in the project as git submodules. The easiest way to initialise these is with the `easy_build.sh` script. However, if you are already have some of them available you may wish to not run the script and instead configure them yourself. This is discussed further below. The libraries are:

- [Args](http://github.com/taywee/args)

    An excellent header-only command-line argument parser.

- [Cereal](https://github.com/USCiLab/cereal)

    A serialization library for C++. Can write to multiple formats (e.g. XML) but QUIT exclusively uses JSON to exploit a couple of useful features.

- [Eigen](http://eigen.tuxfamily.org)

    A high performace linear algebra library. Used for all signal equations in QUIT.

- [Ceres](http://ceres-solver.org)

    A high performace, high quality collection of non-linear least squares solvers.

- [ITK](http://itk.org)

    The Insight Tool-Kit. This is intended as a registration library for medical images. Although the registration features are not used in QUIT, the overall framework is extremely useful. In particular, ITK reads a large variety of common file-formats, and all ITK filters automatically check that their inputs share a common space / co-ordinate definition. This single feature alone prevents the vast majority of easy user mistakes - e.g. trying to use a B1 map that has a different voxel spacing to the input data.

## Compilation

If you are unfamiliar with C++/CMake/git etc. a script is provided that should be able to build the tools provided the right software is available on your system. To use it, in a terminal window, change directory to where you unpacked QUIT. Then type `./easy_build.sh`. This should correctly checkout the git repositories for each external library, build Ceres and ITK, and then configure and build QUIT.

If you are familiar with C++/CMake/git etc. then compiling QUIT should be straightforward. The rough steps are:

1. Ensure you have all the libraries above available. If you do not have system versions, then run `git submodule init; git submodule update` in the QUIT directory.
2. Compile the Ceres library. Make sure it uses the same version of Eigen that you will use for QUIT.
3. Compile ITK following their instructions.
4. Create a build directory for QUIT and use cmake/ccmake to configure the project. Specify the paths to the library when prompted. It is likely that CMake will report an error on the first 'Configure' attempt if it cannot locate Eigen, specify the correct directory and run the configure step again.
5. Compile QUIT.

## File Formats

The available file formats are controlled by the main `CMakeLists.txt` in the root QUIT directory, by listing them as `COMPONENTS` in the ITK `find_package()` step. Add any additional file formats you wish to use here.s

## Tests

QUIT programs are tested using the Bash Automated Test System [(BATS)](http://github.com/bats-core/bats-core), which is included as a git submodule. To run the tests, point BATS at the Tests directory, e.g. `path/to/bats QUIT/Test`, which will run all of the tests. It is possible to run the test files individually as well. It was decided to use BATS instead of a unit-test based framework because this allows the QUIT programs to be tested as a whole, including command-line arguments.

Most QUIT programs are tested by generating ground-truth parameter files with `qinewimage`, feeding these into `qisignal` to generate simulated MR images with added noise, and then running the particular QUIT program to calculate some parameter maps, then comparing these to the ground-truth with `qidiff`. `qidiff` calculates figure-of-merit based on noise factors, i.e. they are a measure of how much the signal noise is amplified in the final maps. In this way the tests also serve to illustrate the quality of the methods as well as whether the programs run correctly. For programs where a ground-truth image cannot be generated easily, the tests at least ensure that the program runs and does not crash.

## The ApplyAlgorithmFilter

The core part of QUIT is the `ApplyAlgorithmFilter` and its child-class `Algorithm`, found in `/Source/Filters/`. This is a sub-class of the ITK `ImageToImageFilter`. The vast majority of QUIT programs declare an `Algorithm` sub-class and use this to process the data. `ApplyAlgorithmFilter` abstracts out most of the heavy lifting of extracting voxel-wise data from multiple inputs and writing it out to multiple outputs, leaving the `Algorithm` to process a single-voxel.

An `Algorithm` defines the number of expected inputs and their size, the number of 'constants' or fixed-parameters, and the number of outputs. It would be preferable if `Algorithm` also defined the types of these, and then `Algorithm` was passed as a template-type to `ApplyAlgorithmFilter`, e.g. `ApplyAlgorithmFilter<DESPOT1Algorithm>`. However, due to an `itk::Image<itk::VariableLengthVector, 3>` being different to an `itk::VectorImage<float, 3>` this is not possible. Instead, `ApplyAlgorithmFilter` takes the input and output types as template parameters, and defines a child-class that has these types available to it. It is these child-classes that developers should sub-class. Several are predefined in the `ApplyTypes.h` file.

## Example: qidespot1

The structure of `qidespot1` is similar to most QUIT programs, and is a good example of most features. At the start are the includes (obviously). After that several `Algorithm` subclasses are defined, as well as a Ceres cost-function. The Ceres documentation is excellent, so refer to that for more information. After all the `Algorithm` classes are defined, the main program body begins. At the start of the program, all the command-line options are defined and then parsed. Then the various inputs are read and passed to the `ApplyAlgorithmFilter`, which is then updated. Finally, the outputs are written back to disk.
