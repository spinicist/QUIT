# Developer

Below are a few introductory notes for anyone who might be interested in contributing to QUIT. Although I hope that most QUIT code is fairly clear, the following is provided to give a high-level overview of how most QUIT programs are structured.

QUIT is built using several C++ libraries. These are:

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

## The ApplyAlgorithmFilter

The core part of QUIT is the `ApplyAlgorithmFilter` and its child-class `Algorithm`, found in `/Source/Filters/`. This is a sub-class of the ITK `ImageToImageFilter`. The vast majority of QUIT programs declare an `Algorithm` sub-class and use this to process the data. `ApplyAlgorithmFilter` abstracts out most of the heavy lifting of extracting voxel-wise data from multiple inputs and writing it out to multiple outputs, leaving the `Algorithm` to process a single-voxel.

An `Algorithm` defines the number of expected inputs and their size, the number of 'constants' or fixed-parameters, and the number of outputs. It would be preferable if `Algorithm` also defined the types of these, and then `Algorithm` was passed as a template-type to `ApplyAlgorithmFilter`, e.g. `ApplyAlgorithmFilter<DESPOT1Algorithm>`. However, due to an `itk::Image<itk::VariableLengthVector, 3>` being different to an `itk::VectorImage<float, 3>` this is not possible. Instead, `ApplyAlgorithmFilter` takes the input and output types as template parameters, and defines a child-class that has these types available to it. It is these child-classes that developers should sub-class. Several are predefined in the `ApplyTypes.h` file.

## Example: qidespot1

The structure of `qidespot1` is similar to most QUIT programs, and is a good example of most features. At the start are the includes (obviously). After that several `Algorithm` subclasses are defined, as well as a Ceres cost-function. The Ceres documentation is excellent, so refer to that for more information. After all the `Algorithm` classes are defined, the main program body begins. At the start of the program, all the command-line options are defined and then parsed. Then the various inputs are read and passed to the `ApplyAlgorithmFilter`, which is then updated. Finally, the outputs are written back to disk.