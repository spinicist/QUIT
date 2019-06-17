Developer
=========

Below are a few introductory notes for anyone who wants to compile QUIT from source or might be interested in contributing to QUIT. Although I hope that most QUIT code is fairly clear, the following is provided to give a high-level overview of how most QUIT programs are structured.

Basic Requirements
------------------

1. A C++17 compliant compiler (GCC 7.0.0+, Clang 3.9+)
2. CMake version 3.10.2 or higher (http://www.cmake.org)

WARNING - You will require a recent compiler for QUIT as it uses C++17 features. On Mac simply having the most recent software updates is enough, on Linux you will need GCC 7.0.0. No one has been brave enough to compile QUIT on Windows to date.

Installing GCC 7
----------------

Most Linux systems currently ship with GCC 4.8 or lower as their system compiler. As of QUIT 2.1 you will require GCC 7 as QUIT uses C++17 features. Joost Kuijer has kindly provided the following recipe for installing newer GCC versions in a user directory and using it to compile QUIT.

1. Follow the GCC guide here: https://gcc.gnu.org/wiki/InstallingGCC
2. Before running the ``build.sh`` script let cmake know about the new compiler
    .. code-block:: bash

        export PATH=$HOME/GCC-7.0.0/bin:$PATH
        export LD_LIBRARY_PATH=$HOME/GCC-7.0.0/lib:$HOME/GCC-7.0.0/lib64:$LD_LIBRARY_PATH
        export CC=$HOME/GCC-7.0.0/bin/gcc
        export CXX=$HOME/GCC-7.0.0/bin/g++

3. Before executing the compiled code also do (you can add this line to your `.bashrc` file:
    .. code-block:: bash

        export LD_LIBRARY_PATH=$HOME/GCC-7.0.0/lib:$HOME/GCC-7.0.0/lib64:$LD_LIBRARY_PATH

External Libraries
------------------

QUIT is built using several C++ libraries. These are currently included in the project as git submodules. The easiest way to initialise these is with the ``build.sh`` script. However, if you are already have some of them available you may wish to not run the script and instead configure them yourself. This is discussed further below. The libraries are:

- `Eigen <http://eigen.tuxfamily.org>`_

    A high performace linear algebra library. Used for all signal equations in QUIT.

- `Ceres <http://ceres-solver.org>`_

    A high performace, high quality collection of non-linear least squares solvers.

- `ITK <http://itk.org>`_

    The Insight Tool-Kit. This is intended as a registration library for medical images. Although the registration features are not used in QUIT, the overall framework is extremely useful. In particular, ITK reads a large variety of common file-formats, and all ITK filters automatically check that their inputs share a common space / co-ordinate definition. This single feature alone prevents the vast majority of easy user mistakes - e.g. trying to use a B1 map that has a different voxel spacing to the input data.

Compilation
-----------

If you are unfamiliar with C++/CMake/git etc. a script is provided that should be able to build the tools provided the right software is available on your system. To use it, in a terminal window, change directory to where you unpacked QUIT. Then type ``./build.sh``. This should correctly checkout the git repositories for each external library, build Ceres and ITK, and then configure and build QUIT.

CMake projects separate the ``build`` and ``install`` phases. The binaries are only moved to a single folder during ``install``. By default, the install directory is ``/usr/local/bin``. If you run ``./build.sh -i``it will run the install step to this directory. You can also ``cd`` into the ``build/`` directory and then type ``make install`` (or ``ninja install``) to avoid re-running the entire ``build.sh`` script.
 
If you want to change where the binaries are installed, you can either:
- Run ``./build.sh -i -p /path/to/install``
- ``cd /build; ccmake ./``, change the ``INSTALL_PREFIX_DIR``, configure (press c), generate (press g), exit, then ``make install``
 
Note that the binaries will end up in ``/bin`` inside the install prefix.

If you are familiar with C++/CMake/git etc. then compiling QUIT should be straightforward. The rough steps are:

1. Ensure you have all the libraries above available. If you do not have system versions, then run ``git submodule init; git submodule update`` in the QUIT directory.
2. Compile the Ceres library. Make sure it uses the same version of Eigen that you will use for QUIT.
3. Compile ITK following their instructions.
4. Create a build directory for QUIT and use cmake/ccmake to configure the project. Specify the paths to the library when prompted. It is likely that CMake will report an error on the first 'Configure' attempt if it cannot locate Eigen, specify the correct directory and run the configure step again.
5. Compile QUIT.

File Formats
------------

The available file formats are controlled by the main ``CMakeLists.txt`` in the root QUIT directory, by listing them as ``COMPONENTS`` in the ITK ``find_package()`` step. Add any additional file formats you wish to use here.

Tests
-----

QUIT programs are tested using the ``nipype`` wrappers and ``unittest``. To run the tests, first build and install ``QUIT`` (so that the executables are available in your ``$PATH`` variable), then ``cd`` into ``Python/Tests`` and then run ``python -m unittest discover``, which will run all of the tests. It is possible to run the test files individually as well.

Most QUIT programs are tested by generating ground-truth parameter files with ``qinewimage``, feeding these into each ``QUIT`` program with the ``--simulate`` argument to generate simulated MR images with added noise, and then running them back through the ``QUIT`` program to calculate the parameter maps, and comparing these to the ground-truth with ``qidiff``. ``qidiff`` calculates a figure-of-merit based on noise factors, i.e. they are a measure of how much the signal noise is amplified in the final maps. In this way the tests also serve to illustrate the quality of the methods as well as whether the programs run correctly. For programs where a ground-truth image cannot be generated easily, the tests at least ensure that the program runs and does not crash.

The ModelFitFilter
------------------

The core part of QUIT is the ``ModelFitFilter`` and its dependent type ``FitFunction``, found in ``Source/Core/``. This is a sub-class of the ITK ``ImageToImageFilter``. The vast majority of QUIT programs declare an `Model` and `FitFunction` sub-class and use these to process the data. ``ModelFitFilter`` abstracts out most of the heavy lifting of extracting voxel-wise data from multiple inputs and writing it out to multiple outputs, leaving the ``FitFunction`` to process a single-voxel. A ``Model`` defines the number of expected inputs and their size, the number of fixed & varying parameters, and the number of outputs.

Example: ``qidespot1``
----------------------

The structure of ``qidespot1`` is similar to most QUIT programs, and is a good example of most features. At the start are the includes (obviously). After that several ``FitFunction`` subclasses are defined, as well as a Ceres cost-function. The Ceres documentation is excellent, so refer to that for more information. After all the ``FitFunction`` classes are defined, the main program body begins. At the start of the program, all the command-line options are defined and then parsed. Then the various inputs are read and passed to the ``ModelFitFilter``, which is then updated. Finally, the outputs are written back to disk.
