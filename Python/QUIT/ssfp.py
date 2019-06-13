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

from os import path, getcwd
from nipype.interfaces.base import TraitedSpec, File, traits
from . import base as QI

############################ qi_ssfp_bands ############################
# < To be implemented > #

############################ qi_ssfp_elipse ############################


class EllipseInputSpec(QI.FitInputSpec):
    # Additional Options
    algo = traits.String(desc='Choose algorithm (h/d)', argstr='--algo=%s')


class EllipseOutputSpec(TraitedSpec):
    G_map = File('ES_G.nii.gz', desc='Path to G map', usedefault=True)
    a_map = File('ES_a.nii.gz', desc='Path to a map', usedefault=True)
    b_map = File('ES_b.nii.gz', desc='Path to b map', usedefault=True)
    theta0_map = File('ES_theta_0.nii.gz',
                      desc='Path to theta 0 (off-resonance phase) map', usedefault=True)
    phi_rf_map = File('ES_phi_rf.nii.gz',
                      desc='Path to RF phase map', usedefault=True)
    residual_map = File('ES_SoS_residual.nii.gz',
                        desc='Path to residual map', usedefault=True)


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


class eMTInputSpec(QI.InputSpec):
    # Inputs
    G_map = File(exists=True, argstr='%s', mandatory=True,
                 position=-3, desc='Path to G parameter map')
    a_map = File(exists=True, argstr='%s', mandatory=True,
                 position=-2, desc='Path to a parameter map')
    b_map = File(exists=True, argstr='%s', mandatory=True,
                 position=-1, desc='Path to b parameter map')

    # Options
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    f0map_file = File(desc='f0 map (Hertz) file', argstr='--f0=%s')
    T2b = traits.Float(desc='T2 of bound pool', argstr='--T2b=%f')


class eMTOutputSpec(TraitedSpec):
    PD_map = File('EMT_PD.nii.gz', desc="Path to PD map", usedefault=True)
    f_b_map = File('EMT_f_b.nii.gz',
                   desc="Path to bound-pool fraction map", usedefault=True)
    k_bf_map = File('EMT_k_bf.nii.gz',
                    desc="Path to exchange rate map", usedefault=True)
    T1_f_map = File('EMT_T1_f.nii.gz',
                    desc="Path to free-pool T1 map", usedefault=True)
    T2_f_map = File('EMT_T2_f.nii.gz',
                    desc="Path to free-pool T2 map", usedefault=True)
    residual_map = File('EMT_SoS_residual.nii.gz',
                        desc="Path to residual map", usedefault=True)


class eMT(QI.FitCommand):
    """
    Fit MT parameters to ellipse parameters

    """

    _cmd = 'qi_ssfp_emt'
    input_spec = eMTInputSpec
    output_spec = eMTOutputSpec


class eMTSimInputSpec(QI.SimInputBaseSpec):
    # Options
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    f0map_file = File(desc='f0 map (Hertz) file', argstr='--f0=%s')
    T2b = traits.Float(desc='T2 of bound pool', argstr='--T2b=%f')
    G_file = File(argstr='%s', mandatory=True,
                  position=-3, desc='Output G file')
    a_file = File(argstr='%s', mandatory=True,
                  position=-2, desc='Output a file')
    b_file = File(argstr='%s', mandatory=True,
                  position=-1, desc='Output b file')


class eMTSimOutputSpec(TraitedSpec):
    G_file = File(desc='Output G file')
    a_file = File(desc='Output a file')
    b_file = File(desc='Output b file')


class eMTSim(QI.SimCommand):
    """
    Simulate ellipse parameters from MT parameters

    """

    _cmd = 'qi_ssfp_emt'
    _param_files = ['PD', 'f_b', 'k_bf', 'T1_f', 'T2_f']
    input_spec = eMTSimInputSpec
    output_spec = eMTSimOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['G_file'] = path.abspath(self.inputs.G_file)
        outputs['a_file'] = path.abspath(self.inputs.a_file)
        outputs['b_file'] = path.abspath(self.inputs.b_file)
        return outputs

############################ qi_ssfp_planet ############################


class PLANETInputSpec(QI.InputSpec):
    # Inputs
    G_map = File(exists=True, argstr='%s', mandatory=True,
                 position=-3, desc='Path to G parameter map')
    a_map = File(exists=True, argstr='%s', mandatory=True,
                 position=-2, desc='Path to a parameter map')
    b_map = File(exists=True, argstr='%s', mandatory=True,
                 position=-1, desc='Path to b parameter map')

    # Options
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    residuals = traits.Bool(
        desc='Write out residuals for each data-point', argstr='--resids')


class PLANETOutputSpec(TraitedSpec):
    PD_map = File('PLANET_PD.nii.gz', desc="Path to PD map", usedefault=True)
    T1_map = File('PLANET_T1.nii.gz', desc="Path to T1 map", usedefault=True)
    T2_map = File('PLANET_T2.nii.gz', desc="Path to T2 map", usedefault=True)
    residual_map = File('PLANET_SoS_residual.nii.gz',
                        desc="Path to residual map", usedefault=True)


class PLANET(QI.FitCommand):
    """
    Calculate T1/T2 from ellipse parameters

    """

    _cmd = 'qi_ssfp_planet'
    input_spec = PLANETInputSpec
    output_spec = PLANETOutputSpec


class PLANETSimInputSpec(QI.SimInputBaseSpec):
    G_file = File(argstr='%s', mandatory=True,
                  position=-3, desc='Output G file')
    a_file = File(argstr='%s', mandatory=True,
                  position=-2, desc='Output a file')
    b_file = File(argstr='%s', mandatory=True,
                  position=-1, desc='Output b file')
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')


class PLANETSimOutputSpec(TraitedSpec):
    G_file = File(desc='Output G file')
    a_file = File(desc='Output a file')
    b_file = File(desc='Output b file')


class PLANETSim(QI.SimCommand):
    """
    Simulate ellipse parameters from T1/T2

    """

    _cmd = 'qi_ssfp_planet'
    _param_files = ['PD', 'T1', 'T2']
    input_spec = PLANETSimInputSpec
    output_spec = PLANETSimOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['G_file'] = path.abspath(self.inputs.G_file)
        outputs['a_file'] = path.abspath(self.inputs.a_file)
        outputs['b_file'] = path.abspath(self.inputs.b_file)
        return outputs
