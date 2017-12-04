---
title: 'QUIT: QUantitative Imaging Tools"
tags:
  - neuroimaging
  - mri
  - cpp
  - quantitative imaging
authors:
  - Tobias C Wood
    orcid: 0000-0001-7640-5520
    affiliation: 1
affiliations:
  - name: Department of Neuroimaging, King's College London
  - index: 1
date: 2017-11-23
bibliography: paper.bib
---

# Summary

- Magnetic Resonance Imaging (MRI) is evolving from acquiring simple qualitative
images towards quantitative data, where information can be inferred about tissue
within each voxel individually.
- Many papers written in this area do not make processing code available. This
hampers adoption of valuable quantitative techniques outside the groups where
they are developed. There is hence a need for open-source implementations of
common quantitative methods.
- The recent qMRLab project (https://github.com/neuropoly/qMRLab) provides a
large variety of techniques written in Matlab. Although Matlab is widespread
in academic environments, it is not always available. In addition, qMRLab is
unsuited to high-throughput processing of multiple datasets.
- QUIT (http://github.com/spinicist/QUIT) provides implementations of multiple 
quantitative MRI algorithms as command-line tools written in C++ and
using the high-quality Insight Toolkit (http://itk.org) and Eigen 
(http://eigen.tuxfamily.org) libraries. They are multi-threaded where-ever
possible. They are easy to include into shell-script based processing pipelines,
and are easy to use with a queue system such as the Sun Grid Engine. They are
hence suitable for the high-throughput processing that is required for timely
analysis of a large neuroimaging study.

- The original focus of QUIT was T1 & T2 mapping using the DESPOT family of
techniques. It has since expanded to include B1 mapping (AFI & DREAM), standard
multi-echo T2 measurements, Laplacian phase-unwrapping and basic Chemical
Exchange Saturation Transfer (CEST) processing. In addition, it includes a
number of useful utilities, such as Region-Of-Interest value extraction and
generation of design matrices suitable for use with FSL Randomise.
- All QUIT programs will give usage instructions when invoked with the `--help`
argument, or with no arguments. Example scripts that include pre-processing
using FSL are included in the repository, and a Wiki is provided with pages for
the most used programs.
- Because QUIT uses ITK, the structure of the programs should be familiar to
anyone acquainted with the library. ITK is principally oriented at image
registration, which is a related but separate topic to quantitative image
processing. However, the wide variety of input file-formats supported by ITK
and the well thought-out processing structure were significant advantages during
development of QUIT. Most QUIT programs are structured around a custom ITK
filter, the `ApplyAlgorithmFilter`, which abstracts out the process of loading
and processing multiple images on a voxel-wise basis. Hence all that is required
in most of the QUIT programs is writing a small `Algorithm` class which is used
to process each voxel. This structure requires a small amount of boiler-plate
code to be written, but otherwise massively simplifies the code that has to be
implemented for each method.

- QUIT has been used during the preparation of the following manuscripts:
[@Wood:2016b]
[@Hawkins:2017]
[@Richetto:2017]
and is currently in use for several ongoing studies.
