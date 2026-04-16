Perfusion
=========

Perfusion is the study of blood flow within the brain. This module contains tools to calculate Cerebral Blood Flow (CBF), the Oxygen Extraction Fraction (OEF), and Z-shimming.

The following commands are available:

* `qi asl`_
* `qi ase_oef`_
* `qi zshim`_

qi asl
------

This command implements the standard equation to calculate Cerebral Blood Flow (CBF) in units of ml/100 g/minute from either Continuous or pseudo-Continuous Arterial Spin Labelling data (CASL or pCASL). For the exact equation used, see the first reference below.

.. image:: cbf.png
    :alt: A Rodent CBF Map

**Example Command Line**

.. code-block:: bash

    qi asl asl_file.nii.gz --blood=2.429 --alpha=0.9 --average --slicetime --pd=reference_file.nii.gz <input.json

The input file must contain pairs of label & control volumes. Currently the order of these is hard-coded to label, then control. The file can contain multiple pairs if you are studying timeseries data. The arguments are discussed further below. It is highly recommended to provide either a separated Proton Density reference image or a tissue T1 map.

**Example JSON File**

.. code-block:: json

    {
        "CASL" : {
            "TR" : 4.0,
            "label_time" : 3.0,
            "post_label_delay" : [ 0.3, 0.6, 0.9 ]
        }
    }


The units for all these values must be consistent, seconds are preferred. If single-slice or 3D data was acquired, then ``post_label_delay`` should contain a single value. For multi-slice data, specify the ``--slicetime`` option and then provide the effective post-labelling delay for each slice.

**Outputs**

* ``CASL_CBF`` - The CBF value, given in mL/(100 g)/min

*Important Options*

* ``--blood, -b``

    The T1 value of blood at the field strength used. The default value is 1.65 seconds, corresponding to 3T. For 1.5T the value should be 1.35 seconds and at 9.4T it should 2.429 seconds. See reference 2.

* ``--pd, -p``

    Provide a separate image to estimate of the Proton Density of tissue. If this is not provided, the label images are used instead.

* ``--tissue, -t``

    Provide a T1 map to correct the Proton Density estimate. If a separate PD reference is not given, then an alternative is to correct the label images for incomplete T1 relaxation.

* ``--alpha, -a``

    The labelling efficiency of the sequence. Default is 0.9.

* ``--lambda, -l``

    The blood-brain partition co-efficient, default 0.9 mL/g.

* ``--average, -a``

    Average the time-series to produce a single CBF output.

* ``--slicetime, -s``

    Apply slice-time correction. The number of post-label delays must match the number of slices.

* ``--dummies, -d``

    Discard this many dummy pairs from the start of the input. Default is 0.

* ``--mask, -m``

    Only process voxels within the specified mask.

**References**

- `ISMRM Consortium Recommendations <http://dx.doi.org/10.1002/mrm.25197>`_
- `High-field blood T1 times <http://dx.doi.org/10.1016/j.mri.2006.10.020>`_

qi ase_oef
----------

Estimates the Oxygen Extraction Fraction (OEF) from Asymmetric Spin-Echo (ASE) data using the Blockley model. If the signal evolution each side of a spin-echo in the presence of blood vessels is observed carefully, it does not display simple monoexponential T2* decay close to the echo, but is instead quadratically exponential. By measuring the T2* decay in the linear regime using an ASE sequence, it is possible to extrapolate back to the echo and obtain an estimate of what the signal would be if no blood was presence. The difference between this and the observed signal can be attributed to the Deoxygenated Blood Volume (DBV), and from there the OEF can be calculated.

**Example Command Line**

.. code-block:: bash

    qi ase_oef ase_file.nii.gz --B0=9.4 <input.json

**Example JSON File**

.. code-block:: json

    {
        "MultiEcho" : {
            "TR" : 2.0,
            "TE1" : 0,
            "ESP" : 0.002,
            "ETL" : 10
        }
    }

Or:

.. code-block:: json
    {
        "MultiEcho" : {
            "TR" : 2.0,
            "TE" : [-0.004, -0.002, 0.0, 0.002, 0.004, 0.006]
        }
    }

``TR`` must be provided but is not used in the calculation. All echo-times are used in the fit.

**Outputs**

* ``ASE_S0.nii.gz`` - Scaling factor.
* ``ASE_dT.nii.gz`` - Echo time offset.
* ``ASE_R2p.nii.gz`` - The R2 prime map. Units are the same as those used for ``TR``, ``TE1`` and ``ESP``.
* ``ASE_DBV.nii.gz`` - The Deoxygenated Blood Volume (only when ``--DBV`` is not specified).
* ``ASE_Tc.nii.gz`` - The critical time (derived parameter).
* ``ASE_OEF.nii.gz`` - The Oxygen Extraction Fraction (derived parameter).
* ``ASE_dHb.nii.gz`` - The Deoxyhaemoglobin concentration (derived parameter).

*Important Options*

* ``--B0, -B``

    Field-strength the data was acquired at. Default is 3.0 (3 Tesla). This is used to calculate Tc and appears elsewhere in several equations.

* ``--Hct, -h``

    Hematocrit value. Default is 0.34.

* ``--DBV, -d``

    Fix the DBV to the specified value and only fit R2'. When specified, the ``ASE_DBV`` output is not written. Default is 0.0 (DBV is fitted as a free parameter).

**References**

- `Blockley <https://doi.org/10.1016/j.neuroimage.2016.11.057>`_

qi zshim
--------

Performs Z-shimming to compensate for signal loss due to through-slice dephasing. Combines data acquired with multiple Z-shim gradients using a root-sum-of-squares approach.

**Example Command Line**

.. code-block:: bash

    qi zshim zshimmed_file.nii.gz --zshims=8 --yshims=1

Does not read input from ``stdin``. The input file should contain data acquired with multiple Z-shim (and optionally Y-shim) gradients, concatenated into a single 4D file.

**Outputs**

* ``{input}_zshim.nii.gz`` - The combined image after Z-shimming.

*Important Options*

* ``--zshims, -z``

    Number of Z-shims. Default is 8.

* ``--yshims, -y``

    Number of Y-shims. Default is 1.

* ``--zdrop``

    Number of Z-shims to drop from each end. Default is 0.

* ``--ydrop``

    Number of Y-shims to drop from each end. Default is 0.

* ``--noiseregion, -n``

    Specify a region to measure and subtract noise. The argument is a region specification string.

* ``--threads, -T``

    Use N threads (default is hardware limit or ``$QUIT_THREADS``).

* ``--out, -o``

    Add a prefix to output filenames.
