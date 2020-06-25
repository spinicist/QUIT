#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT MT modules.

Contains wrappers for:
    - qi_asl
    - qi_ase_oef
    - qi_zshim

Requires that the QUIT tools are in your your system path

"""

from os import path
from nipype.interfaces.base import CommandLine, TraitedSpec, File, traits, isdefined
from .. import base

############################ qi_asl ############################


class ASLInputSpec(base.FitInputSpec):
    # Additional Options
    average = traits.Bool(
        desc='Average across time-series', argstr='--average')
    dummies = traits.Int(
        desc='Remove dummies from timeseries', argstr='--dummies=%d')
    T1_blood = traits.Float(
        desc='Value of blood T1 (default 1.65)', argstr='--blood=%f')
    T1_tissue = traits.String(
        desc='Path to tissue T1 map (units are seconds)', argstr='--tissue=%s')
    PD_map = traits.String(desc='Path to PD weighted image', argstr='--pd=%s')
    label_efficiency = traits.Float(
        desc='Labelling efficiency, default 0.9', argstr='--alpha=%f')
    blood_brain_partition = traits.Float(
        desc='Blood-Brain Partition Co-efficient, default 0.9 mL/g', argstr='--lambda=%f')
    slicetime = traits.Bool(
        desc='Apply slice-timing correction (number of PLDs and slices must match)', argstr='--slicetime')


class ASL(base.FitCommand):
    """
    Calculate CBF from CASL data

    """

    _cmd = 'qi asl'
    input_spec = ASLInputSpec
    output_spec = base.FitOutputSpec('CASL', 'CBF')

############################ qi_ase_oef ############################


class ASEInputSpec(base.FitInputSpec):
    # Additional Options
    B0 = traits.Float(
        desc='Field-strength (Tesla), default 3', argstr='--B0=%f')
    fix_DBV = traits.Float(
        desc='Fix Deoxygenated Blood Volume to value (fraction)', argstr='--DBV=%f')
    gradz = traits.String(
        desc='Field gradient in z/slice-direction for MFG correction', argstr='--gradz=%s')
    slice_thickness = traits.Float(
        desc='Actual slice-thickness for MFG correction if slice-gap was used', argstr='--slice=%f')


class ASEOutputSpec(TraitedSpec):
    S0_map = File('ASE_S0.nii.gz', desc='Path to S0', usedefault=True)
    dT_map = File('ASE_dT.nii.gz',
                  desc='Path to offset time map (s)', usedefault=True)
    R2p_map = File('ASE_R2p.nii.gz',
                   desc='Path to R2′ map (s)', usedefault=True)
    DBV_map = File(
        'ASE_DBV.nii.gz', desc='Path to Deoxygenated Blood Volume map (fraction)', usedefault=True)
    Tc_map = File(
        'ASE_Tc.nii.gz', desc='Path to characteristic time (Tc) map (s)', usedefault=True)
    OEF_map = File(
        'ASE_OEF.nii.gz', desc='Path to Oxygen Extraction Fraction map (fraction)', usedefault=True)
    dHb_map = File(
        'ASE_dHb.nii.gz', desc='Path to Deoxyhaemoglobin concentration map', usedefault=True)


class ASE(base.FitCommand):
    """
    Calculate the Oxygen Extraction Fraction from Asymmetric Spin Echo data

    """

    _cmd = 'qi ase_oef'
    input_spec = ASEInputSpec
    output_spec = base.FitOutputSpec(
        'ASE', ['S0', 'dT', 'R2p', 'DBV', 'Tc', 'OEF', 'dHb'])


class ASESim(base.SimCommand):
    """
    Simulate the Oxygen Extraction Fraction from Asymmetric Spin Echo data with fixed DBV
    """

    _cmd = 'qi ase_oef'
    input_spec = base.SimInputSpec('ASE',
                                   varying=['S0', 'dT', 'R2p'],
                                   extras={'B0': traits.Float(desc='Field-strength (Tesla), default 3', argstr='--B0=%f'),
                                           'fix_DBV': traits.Float(desc='Fix Deoxygenated Blood Volume to value (fraction)', argstr='--DBV=%f', mandatory=True)})
    output_spec = base.SimOutputSpec('ASE')


class ASEDBVSim(base.SimCommand):
    """
    Simulate OEF with a varying DBV
    """
    _cmd = 'qi ase_oef'
    input_spec = base.SimInputSpec('ASE',
                                   varying=['S0', 'dT', 'R2p', 'DBV'],
                                   extras={'B0': traits.Float(desc='Field-strength (Tesla), default 3', argstr='--B0=%f')})
    output_spec = base.SimOutputSpec('ASE')

############################ qi_zshim ############################


class ZShimInputSpec(base.InputBaseSpec):
    in_file = File(argstr='%s', mandatory=True, exists=True,
                   position=-1, desc='Input file to fit polynomial to')
    zshims = traits.Int(argstr='--zshims=%d', desc='Number of Z-shims')
    yshims = traits.Int(argstr='--yshims=%d', desc='Number of Y-shums')
    prefix = traits.String(
        argstr='--out=%s', desc='Output prefix (default is input filename)')
    noiseregion = traits.String(
        desc='Noise region for background correction', argstr='--noiseregion=%s')


class ZShimOutputSpec(TraitedSpec):
    out_file = File(desc="Shimmed Image")


class ZShim(base.BaseCommand):
    """
    Combine an EPI image with Z/Y-shimming
    """
    _cmd = 'qi zshim'
    input_spec = ZShimInputSpec
    output_spec = ZShimOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        if isdefined(self.inputs.prefix):
            fname = self._gen_fname(self.inputs.prefix, suffix='_zshim')
        else:
            fname = self._gen_fname(self.inputs.in_file, suffix='_zshim')
        outputs['out_file'] = fname
        return outputs
