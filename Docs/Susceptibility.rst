Susceptibility
==============

Susceptibility is a fundamental magnetic property of a material, and determines whether materials are paramagnetic (positive susceptibility) or diamagnetic (negative susceptibility). Quantitative Susceptibility Mapping (QSM) is a branch of MRI that aims to measure the susceptiblity of objects from the phase of the MR data. QUIT currently does not contain a full QSM processing pipeline, but does contain some phase unwrapping tools.

* `qi_unwrap_path`_
* `qi_unwrap_laplace`_

qi_unwrap_path
--------------

An implementation of the quality-guided path-based unwrapping of Abdul-Rahman et al. This is the recommended method to use (preferable over Laplacian).

**Example Command Line**

.. code-block:: bash

    qi_unwrap_path phase_file.nii.gz

The phase file must be specified in radians (i.e. between :math:`-\pi` and :math:`+\pi`). Does not read input from ``stdin``, and currently there are no arguments to control the algorithms behaviour.

**Outputs**

* ``input_unwrapped.nii.gz`` - The unwrapped phase value, in radians.

**References**

- `Abdul-Rahman et al 1 <http://ao.osa.org/abstract.cfm?URI=ao-46-26-6623>`_
- `Abdul-Rahman et al 2 <http://ao.osa.org/abstract.cfm?URI=ao-48-23-4582>`_

qi_unwrap_laplace
-----------------

Implements Laplacian-based phase-unwrapping. Along with phase-unwrapping, the Laplacian method implicitly removes background fields. This means it can alter phase values in undesirable ways and hence is not the preferred method.

**Example Command Line**

.. code-block:: bash

    qi_unwrap_laplace phase_file.nii.gz

The phase file must be specified in radians (i.e. between -pi and +pi). Does not read input from `stdin`.

**Outputs**

* ``input_unwrapped.nii.gz`` The unwrapped phase, in radians.

**Important Options**

* ``--mask, -m``

    Specify a mask, the phase will only be unwrapped inside this.

* ``--erode, -e``

    Radius to erode the input mask by (default 1 mm).

**References**

- `Bakker et al <http://linkinghub.elsevier.com/retrieve/pii/S0730725X12000124>`_
