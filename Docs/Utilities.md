# Utilities

QUIT contains a number of utilities. Note that these are actually compiled in two separate modules - `CoreProgs` contains the bare minimum of programs for the QUIT tests to run, while the actual `Utils` modules contains a larger number of useful tools for neuro-imaging pipelines. Their documentation is combined here.

* [qi_coil_combine](#qi_coil_combine)
* [qi_rfprofile](#qi_rfprofile)
* [qiaffine](#qiaffine)
* [qicomplex](#qicomplex)
* [qihdr](#qihdr)
* [qikfilter](#qikfilter)
* [qimask](#qimask)
* [qipolyfit](#qipolyfit)
* [qipolyimg](#qipolyimg)
* [qireorder](#qireorder)
* [qisplitsubjects](#qisplitsubjects)
* [qidiff](#qidiff)
* [qinewimage](#qinewimage)
* [qisignal](#qisignal)

## qi_coil_combine

The program implements both the COMPOSER and Hammond methods for coil combination. For COMPOSER, a wrapper script that includes registration and resampling of low resolution reference data to the image data can be found in `qi_composer.sh`.

**Example Command Line**

```bash
qi_coil_combine multicoil_data.nii.gz --composer=composer_reference.nii.gz
```

Both the input multi-coil file and the reference file must be complex valued. Does not read input from `stdin`. If a COMPOSER reference file is not specifed, then the Hammond coil combination method is used.

**Outputs**

* `input_combined.nii.gz` - The combined complex-valued image.

**Important Options**

* `--composer, -c`

    Use the COMPOSER method. The reference file should be from a short-echo time reference scan, e.g. UTE or ZTE. If

* `--coils, -C`

    If your input data is a timeseries consisting of multiple volumes, then use this option to specify the number of coils used in the acquisition. Must match the number of volumes in the reference image. Does not currently work with the Hammond method.


* `--region, -r`

    The reference region for the Hammond method. Default is an 8x8x8 cube in the center of the acquisition volume.

**References**

- [COMPOSER][1]
- [Hammond Method][2]

[1]: http://doi.wiley.com/10.1002/mrm.26093
[2]: http://linkinghub.elsevier.com/retrieve/pii/S1053811907009998

## qi_rfprofile

This utility takes a B1+ (transmit field inhomogeneity) map, and reads an excitation slab profile from `stdin`. The two are multiplied together along the slab direction (assumed to be Z), to produce a relative flip-angle or B1 map.

**Example Command Line**

```bash
qihdr b1plus_map.nii.gz output_b1_map.nii.gz < input.txt
```

**Example Input File**

```json
{
    "rf_pos" : [ -5, 0, 5],
    "rf_vals" : [[0, 1, 0],
                 [0, 2, 0]]
}
```

`rf_pos` specifies the positions that values of the RF slab have been calculated at, which are specified in `rf_vals`. Note that `rf_vals` is an array of arrays - this allows `qi_rfprofile` to calculate profiles for multiple flip-angles in a single pass. The units for `rf_pos` are the same as image spacing in the header (usually mm). `rf_vals` is a unitless fraction, relative to the nominal flip-angle.

These values should be generated with a Bloch simulation. Internally, they are used to create a spline to represent the slab profile. This is then interpolated to each voxel's Z position, and the value multiplied by the input B1+ value at that voxel to produce the output.

**Outputs**

* output_b1map.nii.gz - The relative flip-angle/B1 map

## qihdr

Prints the header of input files as seen by ITK to `stdout`. Can extract single header fields or print the entirety.

**Example Command Line**

```bash
qihdr input_file1.nii.gz input_file2.nii.gz --verbose
```

Multiple files can be queried at the same time. The `--verbose` flag will make sure you can tell which is which.

**Important Options**

If any of the following options are specified, then only those fields will be printed instead of the full header. This is useful if you want to use a header field in a script:
* `--origin, -o`
* `--spacing, -S` - The voxel spacing
* `--size, -s` - The matrix size
* `--voxvol, -v` - The volume of one voxel

Another useful option is `--meta, -m`. This will let you query specific image meta-data from the header. You must know the exact name of the meta-data field you wish to obtain.

## qicomplex

Manipulate complex/real/imaginary/magnitude/phase data.

**Example Command Line**

```bash
qicomplex -m input_magnitude.nii.gz -p input_phase.nii.gz -R output_real.nii.gz -I output_imaginary.nii.gz
```

Lower case arguments `--mag, -m, --pha, -p, --real, -r, --imag, -i, --complex, -x` are inputs (of which it is only valid to specify certain combinations, complex OR magnitude/phase OR real/imaginary).

Upper case arguments `--MAG, -M, --PHA, -P, --REAL, -R, --IMAG, -I, --COMPLEX, -X` are outputs, any or all of which can be specified.

An additional input argument, `--realimag` is for Bruker "complex" data, which consists of all real volumes followed by all imaginary volumes, instead of a true complex datatype.

The `--fixge` argument fixes the lack of an FFT shift in the slab direction on GE data by multiplying alternate slices by -1. `--negate` multiplies the entire volume by -1. `--double` reads and writes double precision data instead of floats.




