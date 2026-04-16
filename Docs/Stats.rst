Statistics / GLM Tools
======================

QUIT contains a few tools to help prepare your data for statistical analysis with outside tools, for instance non-parametric tests with `Randomise <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Randomise>`_ or an ROI analysis using `Pandas <https://pandas.pydata.org>`_. These tools are:

* `qi glm_setup`_
* `qi glm_contrasts`_
* `qi rois`_

qi glm_setup
-----------

FSL randomise takes a single 4D file with one volume per subject/timepoint as input, along with some simple text files that represent the GLM. Creating these files can be tedious, particularly with the FSL GUI. This tool makes it quick to create the relevant files.

**Example Command Line**

.. code-block:: bash

    qi glm_setup --groups=groups.txt --covars="brain_volume.txt,brain_volume.txt" --design=glm.txt --out=merged.nii --sort subject_dirs*/D1_T1.nii

This command line will merge all the T1 maps in the directories matching the pattern. The file ``groups.txt`` should contain a single number per line, one for each T1 map. The number represents the group or cell that image belongs to. A group of 0 means exclude this file (so you don't have to work out a pattern that won't match that file). For example, with 8 scans belonding to 3 groups with 1 excluded scan, the ``groups.txt`` might look like:

.. code-block:: bash

    1
    3
    3
    2
    0
    1
    2
    1

The design matrix corresponding to the specified groups will be saved to the ``glm.txt`` file (Note - this will still need to be processed with ``Text2Vest`` to make it compatible with ``randomise``). If ``--sort`` is specified, then the images and design matrix will be sorted into ascending order.

*Important Options*

* ``--groups, -g``

    Path to a file with one group number per line (one per input file). Required. A group of 0 means exclude that file.

* ``--out, -o``

    Path for the output merged 4D file. Required.

* ``--covars, -C``

    Path(s) to covariate files, separated by commas. Each file should have one value per line per input image.

* ``--design, -d``

    Path to save the design matrix.

* ``--contrasts, -c``

    Generate and save auto-generated contrast rows to the specified path. For each group, a contrast row comparing group g to group (g+1) is generated. Covariate contrasts are also generated.

* ``--ftests, -f``

    Generate and save F-test rows to the specified path. One F-test per group, plus a main-effect F-test if there are more than 2 groups.

* ``--sort, -s``

    Sort merged file and design matrix in ascending group order.

qi glm_contrasts
---------------

Randomise does not save any ``contrast`` files, i.e. group difference maps, it only saves the statistical maps. For quantitative imaging, the contrasts can be informative to look at, as if scaled correctly, they can be interpreted as effect size maps. A group difference in human white matter T1 of only tens of milliseconds, even if it has a high p-value, is perhaps not terribly interesting as it corresponds to a change of about 1%. These contrast maps are particularly useful if used with the `dual-coding <https://github.com/spinicist/nanslice>`_ visualisation technique.

**Example Command Line**

.. code-block:: bash

    qi glm_contrasts merged_images.nii design.txt contrasts.txt --out=contrast_prefix

The design and contrasts files should be raw text (not passed through ``Text2Vest``). Comment lines starting with ``#`` or ``/`` are ignored. One contrast image will be generated for each row of the contrast matrix. Output files are named ``{prefix}con1.nii.gz``, ``{prefix}con2.nii.gz``, etc.

*Important Options*

* ``--frac, -F``

    Output contrasts as a fraction of the grand mean.

* ``--out, -o``

    Prefix for output filenames.

qi rois
-------

An alternative to voxel-wise statistics is to average the values over a pre-defined, anatomically meaningful region-of-interest in each quantitative image, and the perform statistics on those ROI values. This approach has several advantages, as more traditional and robust statistical methods can be used than the simple parametric T-tests that voxel-wise analysis tools use.

To avoid resampling issues, it is preferable to warp the ROI definitions (atlas files) to the subject space and sample the quantitative maps at their native resolution. This can make extracting all the ROI values tedious. This tool can extract ROI values from multiple files at once and output a table to ``stdout`` for use with a stats tool such as `Pandas <https://pandas.pydata.org>`_. Redirect output to a file with ``>``. It can also calculate the volumes of the warped ROIs, i.e. for a Tensor/Deformation Based Morphometry analysis. The registrations required for this should be carried out with external tools, e.g. `ANTs <https:://github.com/stnava/ANTs>`_ or `FSL <https://fsl.fmrib.ox.ac.uk>`_.

**Example Command Line**

For quantitative ROIs:

.. code-block:: bash

    qi rois labels_subject1.nii labels_subject2.nii ... labels_subjectN.nii data_subject1.nii data_subject2.nii ... data_subjectN.nii --ignore_zero --header=subject_ids.txt

For ROI volumes:

.. code-block:: bash

    qi rois --volumes labels_subject1.nii labels_subject2.nii ... labels_subjectN.nii --header=subject_ids.txt

Any header files should contain one line per subject, corresponding to the input image files. The output of ``qi rois`` is fairly flexible, and can be controlled with the ``--transpose``, ``--delim`` (default ``,``), ``--precision`` (default 6), and ``--sigma`` options.

*Important Options*

* ``--volumes, -V``

    Output ROI volumes instead of average values. Does not require value images.

* ``--labels, -l``

    Path to a file specifying labels and names. Format is ``label_number,label_name`` one per line.

* ``--print_names, -n``

    Print label names in the first column/row. Requires ``--labels`` to be specified.

* ``--ignore_zero, -z``

    Ignore label 0 (background).

* ``--sigma, -s``

    Print standard deviation along with mean.

* ``--precision, -p``

    Number of decimal places. Default is 6.

* ``--delim, -d``

    Delimiter between entries. Default is ``,``.

* ``--transpose, -t``

    Transpose output table (values in rows instead of columns).

* ``--header, -H``

    Add a header from the specified file (one entry per subject). Can be specified multiple times.

* ``--header_name``

    Header name for column/row labels. Must be specified in the same order as ``--header``.

* ``--scale``

    Divide ROI values by a scale factor. Must match the order of input paths.
