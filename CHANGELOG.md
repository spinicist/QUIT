# QUIT Changelog

## Version 2.0.2

Mostly another bug-fix release, but also to celebrate acceptance by [JOSS](https://doi.org/10.21105/joss.00656)

Notable changes:

1. Added a new Z-shim calculation tool
2. Scaling of residuals by the PD in DESPOT tools has been *removed*

(Yes, this changelog is missing from the actual v2.0.2 tag and downloads. Whoops)

## Version 2.0.1

A bug-fix release

1. JSON-ified the input for `qi_rfprofile`, and added `lower_bounds` and `upper_bounds` to the JSON input for `qimcdespot`
2. Improved the numerical stability of `qipolyfit` so that higher order (8 and above) polynomials fit correctly
3. `qiaffine` did not handle non-isotropic voxel spacings correctly
4. A lot of documentation fixes

## Version 2.0.0

This was a major upgrade from QUIT 1.1 with a lot of improvements and new features. Principal changes were:

1. Introduction of proper [documentation](https://spinicist.github.io/QUIT)
2. A switch to JSON as the input file format
3. Standardisation of input arguments across all tools
4. Introduction of a lot of new methods, including ASL, phase unwrapping and SSFP ellipse methods
