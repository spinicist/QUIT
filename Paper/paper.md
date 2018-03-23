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
date: 2018-03-23
bibliography: paper.bib
---

# Summary

- Magnetic Resonance Imaging (MRI) is evolving from acquiring simple qualitative images towards quantitative data, where information can be inferred about tissue microstructure within each voxel individually.
- Many papers written in this area do not make processing code available. This hampers adoption of valuable quantitative techniques outside the groups where they are developed. There is hence a need for open-source implementations of common quantitative methods.
- QUIT (http://github.com/spinicist/QUIT) provides implementations of multiple quantitative MRI algorithms as command-line tools written in C++ and using the high-quality Insight Toolkit [@ITK], CERES [@CERES] and Eigen [@Eigen] libraries. They are multi-threaded where-ever possible. They are easy to include into shell-script based processing pipelines, and are easy to use with a queue system such as the Sun Grid Engine. They are hence suitable for the high-throughput processing that is required for timely analysis of a large neuroimaging study.
- The recent qMRLab project (https://github.com/neuropoly/qMRLab) has similar aims and provides a large variety of techniques written in Matlab. Although Matlab is widespread in academic environments, it is not always available. In addition, the primary focus of qMRLab is educational and it is unsuited to high-throughput processing of multiple datasets.
- QUIT is currently in use for several ongoing studies at multiple sites and has been used during the preparation of several manuscripts[@Wood:2016b][@Hawkins:2017][@Richetto:2017][@Vanes:2018].

# References
