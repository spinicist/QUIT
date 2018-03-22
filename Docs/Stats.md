# Statistics / GLM Tools

QUIT contains a few tools to help prepare your data for statistical analysis with outside tools, for instance non-parametric tests with [Randomise](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Randomise) or an ROI analysis using [Pandas](https://pandas.pydata.org). These tools are:

* [qi_glmsetup](#qi_glmsetup)
* [qi_glmcontrasts](#qi_glmcontrasts)
* [qi_rois](#qi_rois)

## qi_glmsetup

FSL randomise takes a single 4D file with one volume per subject/timepoint as input, along with some simple text files that represent the GLM. Creating these files can be tedious, particularly with the FSL GUI. This tool makes it quick to create the relevant files.

**Example Command Line**

```bash
qi_glmsetup --groups=groups.txt --covars="brain_volume.txt,brain_volume.txt" --design=glm.txt --out=merged.nii --sort subject_dirs*/D1_T1.nii
```

This command line will merge all the T1 maps in the directories matching the pattern. The file `groups.txt` should contain a single number per line, one for each T1 map. The number represents the group or cell that image belongs to. A group of 0 means exclude this file (so you don't have to work out a pattern that won't match that file). For example, with 8 scans belonding to 3 groups with 1 excluded scan, the `groups.txt` might look like:

```bash
1
3
3
2
0
1
2
1
```

The design matrix corresponding to the specified groups will be saved to the `glm.txt` file (Note - this will still need to be processed with `Text2Vest` to make it compatible with `randomise`). If `--sort` is specified, then the images and design matrix will be sorted into ascending order.

## qi_glmcontrasts

Randomise does not save any `contrast` files, i.e. group difference maps, it only saves the statistical maps. For quantitative imaging, the contrasts can be informative to look at, as if scaled correctly, they can be interpreted as effect size maps. A group difference in human white matter T1 of only tens of milliseconds, even if it has a high p-value, is perhaps not terribly interesting as it corresponds to a change of about 1%. These contrast maps are particularly useful if used with the [dual-coding](https://github.com/spinicist/nanslice) visualisation technique.

**Example Command Line**

```bash
qi_glmcontrasts merged_images.nii design.txt contrasts.txt --out=contrast_prefix
```

The design and contrasts files should be raw text (not passed through `Text2Vest`). One contrast image will be generated for each row of the contrast matrix.

## qi_rois

An alternative to voxel-wise statistics is to average the values over a pre-defined, anatomically meaningful region-of-interest in each quantitative image, and the perform statistics on those ROI values. This approach has several advantages, as more traditional and robust statistical methods can be used than the simple parametric T-tests that voxel-wise analysis tools use.

To avoid resampling issues, it is preferable to warp the ROI definitions (atlas files) to the subject space and sample the quantitative maps at their native resolution. This can make extracting all the ROI values tedious. This tool can extract ROI values from multiple files at once and produce a Comma-Separated Value (.csv) file as output for use with a stats tool such as [Pandas](https://pandas.pydata.org). It can also calculate the volumes of the warped ROIs, i.e. for a Tensor/Deformation Based Morphometry analysis. The registrations required for this should be carried out with external tools, e.g. [ANTs](https:://github.com/stnava/ANTs) or [FSL](https://fsl.fmrib.ox.ac.uk).

**Example Command Line**

For quantitative ROIs:

```bash
qi_rois labels_subject1.nii labels_subject2.nii ... labels_subjectN.nii data_subject1.nii data_subject2.nii ... data_subjectN.nii --ignore_zero --header=subject_ids.txt
```

For ROI volumes:

```bash
qi_rois --volumes labels_subject1.nii labels_subject2.nii ... labels_subjectN.nii ---header=subject_ids.txt
```

Any header files should contain one line per subject, corresponding to the input image files. The output of `qi_rois` is fairly flexible, and can be controlled with the `--transpose`, `--delim`, `--precision`, and `--sigma` options.