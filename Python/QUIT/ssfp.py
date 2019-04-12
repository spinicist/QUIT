#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Implementation of nipype interfaces for QUIT SSFP modules.

To be implemented:
    - qi_ssfp_bands
    - qi_ssfp_elipse
    - qi_ssfp_emt
    - qi_ssfp_planet

Requires that the QUIT tools are in your your system path

"""

from nipype.interfaces.base import TraitedSpec, File, traits
from . import base as QI

############################ qi_ssfp_bands ############################
# < To be implemented > #

############################ qi_ssfp_elipse ############################


class EllipseInputSpec(QI.FitInputSpec):
    # Additional Options
    algo = traits.String(desc="Choose algorithm (h/d)", argstr="--algo=%s")


class EllipseOutputSpec(TraitedSpec):
    G_map = File('ES_G.nii.gz', desc="Path to G map", usedefault=True)
    a_map = File('ES_a.nii.gz', desc="Path to a map", usedefault=True)
    b_map = File('ES_b.nii.gz', desc="Path to b map", usedefault=True)
    residual_map = File('ES_residual.nii.gz',
                        desc="Path to residual map", usedefault=True)


class Ellipse(QI.FitCommand):
    """
    Fit an ellipse to SSFP data

    """

    _cmd = 'qi_ssfp_ellipse'
    input_spec = EllipseInputSpec
    output_spec = EllipseOutputSpec


class EllipseSim(QI.SimCommand):
    """
    Simulate SSFP data from ellipse parameters

    """

    _cmd = 'qi_ssfp_ellipse'
    _param_files = ['G', 'a', 'b', 'theta_0', 'phi_rf']
    input_spec = QI.SimInputSpec
    output_spec = QI.SimOutputSpec

############################ qi_ssfp_emt ############################
# < To be implemented > #

############################ qi_ssfp_planet ############################


class PLANETInputSpec(QI.InputSpec):
    # Inputs
    G_map = File(exists=True, argstr='%s', mandatory=True,
                 position=-3, desc='Path to G parameter map')
    a_map = File(exists=True, argstr='%s', mandatory=True,
                 position=-2, desc='Path to a parameter map')
    b_map = File(exists=True, argstr='%s', mandatory=True,
                 position=-2, desc='Path to b parameter map')

    # Options
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    residuals = traits.Bool(
        desc='Write out residuals for each data-point', argstr='--resids')


class PLANETOutputSpec(TraitedSpec):
    PD_map = File('PLANET_PD.nii.gz', desc="Path to PD map", usedefault=True)
    T1_map = File('PLANET_T1.nii.gz', desc="Path to T1 map", usedefault=True)
    T2_map = File('PLANET_T2.nii.gz', desc="Path to T2 map", usedefault=True)
    residual_map = File('ES_residual.nii.gz',
                        desc="Path to residual map", usedefault=True)


class PLANET(QI.FitCommand):
    """
    Fit an ellipse to SSFP data

    """

    _cmd = 'qi_ssfp_planet'
    input_spec = PLANETInputSpec
    output_spec = PLANETOutputSpec


class PLANETSimInputSpec(QI.SimInputSpec):
    # Options
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')


class PLANETSim(QI.SimCommand):
    """
    Simulate ellipse parameters from T1/T2

    """

    _cmd = 'qi_ssfp_planet'
    _param_files = ['PD', 'T1', 'T2']
    input_spec = PLANETSimInputSpec
    output_spec = QI.SimOutputSpec
