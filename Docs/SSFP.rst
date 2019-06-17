SSFP
====

The Steady-State Free-Precession (SSFP), or more precisely balanced-SSFP (bSSFP), sequence is one of the oldest NMR sequences and can be used to give high SNR MR images with mixed T1/T2 contrast in very short scan time. However, it suffers from banding artefacts in areas of off-resonance which limit its clinical applicability. This module contains a tool for removing those banding artefacts, and then further tools for quantitative mapping using the ellipse signal model.

* `qi_ssfp_bands`_
* `qi_ssfp_ellipse`_
* `qi_ssfp_planet`_
* `qi_ssfp_emt`_


qi_ssfp_ellipse
---------------

The most important result of Xiang & Hoff's Geometric Solution paper was that the SSFP signal equation can be expressed as an ellipse in the complex-plane. Shcherbakova built on this and showed it was possible to recover the ellipse parameters *G*, *a*, *b* from at least six phase-increments. They then proceeded to recover T1 & T2 from the ellipse parameters. This utility calculates the ellipse parameters, and ``qi_ssfp_planet`` then processes those parameters to calculate T1 & T2. A non-linear fit is used instead of the algebraic method used by Shcherbakova et al. This is slower, but robust across all

.. image:: ellipse.png
    :alt: SSFP Ellipse Parameters

**Example Command Line**

.. code-block:: bash

    qi_ssfp_ellipse ssfp_data.nii.gz < input.json

The SSFP file must be complex-valued. At least three pairs of opposing phase-increments are recommended (six images in total).

**Outputs**

- ``ES_G`` - The Geometric Solution point of the ellipse. Influences the overall size of the ellipse. This is called \(M\) in the Hoff and Shcherbakova papers, but it is not a measurable magnetization and hence to distinguish it a different letter is used.
- ``ES_a`` - The ellipse parameter that along with \(G\) controls the ellipse size.
- ``ES_b`` - The ellipse parameter that determines how flat or circular the ellipse is.
- ``ES_theta_0`` - The accrued phase due to off-resonance, divide by :math:`2\pi TE` (or :math:`\pi TR`) to find the off-resonance frequency.
- ``ES_phi_rf`` - The effective phase of the RF pulse.

**References**

- `PLANET <http://dx.doi.org/10.1002/mrm.26717>`_

qi_planet
--------------

Converts the SSFP Ellipse parameters into relaxation times.

**Example Command Line**

.. code-block:: bash

    qi_ssfp_planet ES_G.nii.gz ES_a.nii.gz ES_b.nii.gz

**Outputs**

- ``PLANET_T1.nii.gz`` - Longitudinal relaxation time
- ``PLANET_T2.nii.gz`` - Transverse relaxation time
- ``PLANET_PD.nii.gz`` - Apparent Proton Density

**References**

- `PLANET <http://dx.doi.org/10.1002/mrm.26717>`_

qi_ssfp_emt
-----------

Due to the short TR commonly used with SSFP, at high flip-angles the sequence becomes MT weighted. It is hence possible to extract qMT parameters from SSFP data. More details will be in a forthcoming paper.

**Example Command Line**

.. code-block:: bash

    qi_ssfp_emt ES_G.nii.gz ES_a.nii.gz ES_b.nii.gz

**Outputs**

- ``EMT_T1f.nii.gz`` - Longitudinal relaxation time of the free water bool
- ``EMT_T2f.nii.gz`` - Transverse relaxation time of the free water pool
- ``EMT_M0.nii.gz`` - Apparent Proton Density
- ``EMT_F.nii.gz`` - Bound pool fraction
- ``EMT_kf.nii.gz`` - Forward exchange rate

**References**

- `Bieri et al <http://doi.wiley.com/10.1002/mrm.21056>`_
- `Gloor et al <http://doi.wiley.com/10.1002/mrm.21705>`_
