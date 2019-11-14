Magnetization Transfer
======================

MR voxels often contain complex microstructure with multiple different components or pools, each with unique relaxation properties. It is possible for magnetization to be transferred between these pools via several mechanisms, such as exchange of individual protons or entire molecules, or simple dipolar coupling from molecules that are in close proximity. These mechanisms can be studied in the related fields of Magnetization Transfer (MT) and Chemical Exchange Saturation Transfer (CEST). QUIT contains tools for quantitative MT (`qi qmt`_, `qi ssfp_emt`_), CEST (`qi zspec_interp`_, `qi lorentzian`_) and also for processing MPM data (`qi mtsat`_). Example `nipype workflows <https://github.com/spinicist/QUIT/tree/master/Python/QUIT/workflows/cest.py>`_ for CEST analysis are also provided.

* `qi lineshape`_
* `qi qmt`_
* `qi zspec_interp`_
* `qi lorentzian`_
* `qi dipolar_mtr.sh`_
* `qi ssfp_emt`_
* `qi mtsat`_

qi lineshape
------------

A utility to sample lineshapes and write them out to files, which can then be read by `qi qmt`_ and interpolated values used instead of calculating the lineshape during fitting. This is used principally to speed up fitting the Super-Lorentzian lineshape, which takes a long time to calculate but is smoothly valued, so can be accurately approximated using interpolation.

**Example Command Line**

.. code-block:: bash

    qi lineshape --lineshape=SuperLorentzian --frq_start=500 --frq_space=500 --frq_count=150 > superlorentz.json

No input file required. Output is to ``stdout`` in JSON format, so should be redirected into a file for further use.

**Important Options**

* ``--lineshape, -l``

    Choose from a Gaussian, a Lorentzian or the Super-Lorentzian

* ``--T2b, -t``

    Specify the nominal T2 of the lineshape. During fitting in `qi qmt`_ scaling will be used to find the actual value. Should be specified in seconds.

* ``--frq_count``, ``--frq_start``, ``--frq_space``

    These Control the position and number of samples to take on the lineshape. ``frq_start`` and ``frq_space`` should be in Hertz.

qi qmt
------

Calculates Quantitative Magnetization Transfer parameters using the Ramani model from steady-state gradient-echo data acquired with multiple off-resonance saturation pulses.

.. image:: f_b.png
    :alt: Bound Pool Fraction

**Example Command Line**

.. code-block:: bash

    qi qmt --verbose MTSatData.nii.gz T1.nii.gz --mask brain_mask.nii.gz --lineshape=superlorentz.json --B1=B1_map.nii.gz --f0=B0_map.nii.gz < input.json

Note that a T1 map is a required input to stabilise the fitting.

**Example JSON File**

.. code-block:: json

    {
        "MTSat" : {
            "TR": 0.055,
            "FA": 5,
            "sat_f0": [56360, 47180, 12060, 1000, 1000, 2750, 2770, 2790, 2890, 1000, 1000],
            "sat_angle": [332, 628, 628, 332, 333, 628, 628, 628, 628, 628, 628],
            "pulse": { "name": "Gauss", "Trf": 0.015, "p1": 0.416, "p2": 0.295 }
        }
    }

``Trf`` is the pulse-width (RF Time). ``p1`` and ``p2`` are the ratio of the integral of :math:`B_1` and :math:`B_1^2` (the integrals of the pulse amplitude and the square of the pulse amplitude) to the maximum amplitude of the pulse. Both ``Trf`` and ``TR`` should be in seconds. ``sat_f0`` is in Hertz.

**Outputs**

- ``QMT_f_b.nii.gz`` - The bound pool fraction
- ``QMT_k_bf.nii.gz`` - The forward exchange rate from bound to free pool
- ``QMT_T1_f.nii.gz`` - T1 of the free pool
- ``QMT_T2_f.nii.gz`` - T2 of the free pool
- ``QMT_T2_b.nii.gz`` - T2 of the bound pool
- ``QMT_PD.nii.gz`` - The apparent Proton Density / size of the free pool

Note that ``T1_b``, the longitudinal relaxation rate of the bound pool, is fixed to 1 second in this model and so is not written out.

**References**

- `Ramani et al <http://linkinghub.elsevier.com/retrieve/pii/S0730725X02005982>`_

qi zspec_interp
---------------

Interpolates a Z-spectrum to arbitrary precision. Can output asymmetry values instead of a Z-spectrum.

**Example Command Line**

.. code-block:: bash

    qi zspec_interp zspectrum.nii.gz --f0=LTZ_f0.nii.gz < input.json

The off-resonance map units must match the input frequencies (e.g. either PPM or Hertz)

**Example JSON File**

.. code-block:: json

    {
        "input_freqs" : [ -5, -2.5, 0, 2.5, 5],
        "output_freqs" : [ -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]
    }

``input_freqs`` are the offset frequencies the Z-spectrum was acquired at. ``output_freqs`` are the frequencies you want the asymmetry calculated at.

**Outputs**

* ``{input}_interp.nii.gz`` The interpolated Z-spectrum.

*Important Options*

* ``--f0, -f``

    Specify an off-resonance map. Units must be the same as the input & asymmetry frequencies.

* ``-O, --order``

    The order of Spline interpolation used. Default is 3 (cubic).

* ``-a, --asym``

    Output asymmetry (:math`Z(+f) - Z(-f)`) values.

qi lorentzian
-------------

Fits sums of Lorentzian functions to a Z-spectrum. Highly customisable for the number of desired Lorentzian's and their characteristics.

**Example Command Line**

.. code-block:: bash

    qi lorentzian zspectrum.nii.gz < input.json

The Z-spectrum must be a 4D file with each volume acquired at a different offset frequency.

**Example JSON File**

.. code-block:: json

    {
        'MTSat': {
            'pulse': {
                'p1': 0.4,
                'p2': 0.3,
                'bandwidth': 0.39
            },
            'Trf': 0.02,
            'TR': 4,
            'FA': 5,
            'sat_f0': [0, 1, 2, 3, 4, 5],
            'sat_angle': [180, 180, 180, 180, 180],
        },
        'pools' :
        [
            { 
                'name': 'DS',
                'df0': [0, -2.5, 2.5],
                'fwhm': [1.0, 1.e-6, 3.0],
                'A': [0.2, 1.e-3, 1.0],
                'use_bandwidth': True
            },
            {
                'name': 'MT',
                'df0': [-2.5, -5.0, -0.5],
                'fwhm': [50.0, 35.0, 200.0],
                'A': [0.3, 1.e-3, 1.0]
            }
        ]
    }

The input needs to include both the sequence parameters and the characteristics of the Lorentzian "pools" that you wish to fit. Currently the only important information used from the sequence are the saturation offsets, and optionally the bandwidth of the pulse. For each pool a name is required, and then triples of values representing the starting, lower and upper bound for the center frequency ``df0``, the Full-Width Half-Maximum ``fwhm`` and amplitude ``A`` of the Lorentzian. You can also specify that the modified Lorentzian including the pulse bandwidth should be used `'use_bandwith' : True`. See the reference for details.

*Important Options*

* ``--add, -a``

    Specify an additive model instead of the default subtractive (saturation) model. Useful when a base-line has already been subtracted from the Z-spectrum. See reference for details.

* ``--zref, -z``

    Change the reference value for the Z-spectrum. Default is 1.0, change to 0.0 for additive model.


**Outputs**

For each pool three outputs will be written, prefixed by the pool name. For a single pool representing direct-saturation (DS), the following will be written:

* ``DS_f0.nii.gz``  - The center frequency of the fitted Lorentzian.
* ``DS_fwhm.nii.gz``   - The width of the fitted Lorentzian.
* ``DS_A.nii.gz`` - The amplitdue of the fitted Lorentzian.

**References**

- `Deshmane et al <http://doi.wiley.com/10.1002/mrm.27569>`_

qi dipolar_mtr.sh
-----------------

Calculates dipolar/inhomogeneous Magnetization Transfer Ratios (MTRs). Dipolar/inhomogeneous MT is a new (see note) contrast mechanism that is present in highly structured materials such as myelin and tendon. By applying off-resonance saturation at both positive and negative frequencies (instead of only one side as in classic MTR) it is possible to decouple the dipolar pool and hence produce an enhanced Magnetization Transfer (eMT) effect. The different between eMT and normal MT is the dipolar/inhomogeneous MT and is potentially highly specific to myelin within the brain.

Although the majority of the existing literature refers to this effect as inhomogeneous MT, this name was chosen before the physical phenomena underlying the effect was well understood. Current theory does not rely on inhomogeneous effects at all, so the name is a misnomer.

**Example Command Line**

.. code-block:: bash

    qi dipolar_mtr dipolar_mt_volumes.nii.gz

The input must consist of 5 volumes: Dipolar +/-, Dipolar -/+, Unsaturated, MT+, MT-. This scheme is not flexible and will be improved in a future version.

**Outputs**

* ``DMT_mtr.nii.gz`` - The classic MTR, expressed as a percentage
* ``DMT_emtr.nii.gz`` - The enhanced MTR, expressed as a percentage
* ``DMT_dmtr.nii.gz`` - The dipolar MTR, expressed as a percentage. This is the difference between eMTR and MTR.
* ``DMT_mta.nii.gz`` - The first-order MT-asymmetry (MT- subtracted from MT+, relative to unsaturated, in percent).

**References**

1. `Original full paper <http://doi.wiley.com/10.1002/mrm.25174>`_
2. `Dipolar versus inhomogeneous naming <https://doi.org/10.1016/j.jmr.2016.11.013>`_

qi ssfp_emt
-----------

Due to the short TR commonly used with SSFP, at high flip-angles the sequence becomes MT weighted. It is hence possible to extract qMT parameters from SSFP data. More details will be in a forthcoming paper.

**Example Command Line**

.. code-block:: bash

    qi ssfp_emt ES_G.nii.gz ES_a.nii.gz ES_b.nii.gz

**Outputs**

- ``EMT_T1f.nii.gz`` - Longitudinal relaxation time of the free water bool
- ``EMT_T2f.nii.gz`` - Transverse relaxation time of the free water pool
- ``EMT_M0.nii.gz`` - Apparent Proton Density
- ``EMT_F.nii.gz`` - Bound pool fraction
- ``EMT_kf.nii.gz`` - Forward exchange rate

**References**

- `Bieri et al <http://doi.wiley.com/10.1002/mrm.21056>`_
- `Gloor et al <http://doi.wiley.com/10.1002/mrm.21705>`_

qi mtsat
-----------

Implementation of Gunther Helm's MT-Sat method. Calculates R1, apparent PD and the semi-quantitative MT-Saturation parameter "delta". This is the fractional reduction in the longitudinal magnetization during one TR, expressed as a percentage. Arguably could be included in the :doc:`Relaxometry` module instead. Outputs R1 instead of T1 as this is more common in the MTSat / MPM literature. If using multi-echo input data the input should be passed through `qi mpm_r2s` first and the output ``S0`` files used as input to `qi mtsat`.

**Example Command Line**

.. code-block:: bash

    qi mtsat PDw.nii.gz T1w.nii.gz MTw.nii.gz < input.json

**Example JSON File**

.. code-block:: json

    {
        "MTSat": {
            "TR_PDw": 0.025,
            "TR_T1w": 0.025,
            "TR_MTw": 0.028,
            "FA_PDw": 5,
            "FA_T1w": 25,
            "FA_MTw": 5
        }
    }

**Outputs**

- ``MTSat_R1.nii.gz`` - Apparent longitudinal relaxation rate
- ``MTSat_S0.nii.gz`` - Apparent proton density / equilibrium magnetization
- ``MTSat_delta.nii.gz`` - MT-Sat parameter, see above.

**References**

- `Helms et al <http://doi.wiley.com/10.1002/mrm.21732>`_
- `Erratum <http://doi.wiley.com/10.1002/mrm.22607>`_
