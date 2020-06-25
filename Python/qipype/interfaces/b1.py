#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT B1 utilities
Requires that the QUIT tools are in your your system path
"""

from os import path, getcwd
from nipype.interfaces.base import TraitedSpec, File, traits
from .. import base

############################ qidream ############################
# Implemented but not tested #


class QIDreamInputSpec(base.InputSpec):
    # Inputs
    dream_file = File(exists=True, argstr='%s', mandatory=True,
                      position=0, desc='Input file. Must have 2 volumes (FID and STE)')

    # Options
    threads = traits.Int(
        desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(
        desc='Add a prefix to output filenames', argstr='--out=%s')
    order = traits.String(
        desc='Volume order - f/s/v - fid/ste/vst first', argstr='--order=%s')
    mask_file = File(
        desc='Only process voxels within the mask', argstr='--mask=%s')
    alpha = traits.Float(
        desc="Nominal flip-angle (default 55)", argstr="--alpha=%f")


class QIDreamOutputSpec(TraitedSpec):
    b1_rel_map = File(desc="The relative flip-angle map.")
    b1_act_map = File(desc="The actual achieved angle in each voxel.")


class QIDream(base.BaseCommand):
    """
    Interface for qidream

    Example 1
    -------
    >>> from qipype.interfaces.relax import QiDream
    >>> interface = QiDream(prefix='nipype_', param_file='spgr_params.json')

    """

    _cmd = 'qi dream'
    input_spec = QIDreamInputSpec
    output_spec = QIDreamOutputSpec

    # If the command requires a json input file
    def _format_arg(self, name, spec, value):
        return self._process_params(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['b1_rel_map'] = os.path.abspath(
            self._add_prefix('DREAM_B1.nii.gz'))
        outputs['b1_act_map'] = os.path.abspath(
            self._add_prefix('DREAM_angle.nii.gz'))

        return outputs

############################ qiafi ############################
# Implemented but not tested #


class QIAFIInputSpec(base.InputSpec):
    # Inputs
    afi_file = File(exists=True, argstr='%s', mandatory=True,
                    position=0, desc='Input file')

    # Options
    threads = traits.Int(
        desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(
        desc='Add a prefix to output filenames', argstr='--out=%s')
    alpha = traits.Float(
        desc="Specify nominal flip-angle, default 55", argstr="--flip=%f")
    tr_ratio = traits.Float(
        desc="Specify TR2:TR1 ratio, default 5", argstr="--ratio=%f")
    save_act_b1 = traits.Bool(
        desc="Write out the actual flip-angle as well as B1", argstr="--save")


class QIAFIOutputSpec(TraitedSpec):
    # Specify which outputs there are
    b1_rel_map = File(desc="The relative flip-angle map.")
    b1_act_map = File(desc="The actual flip-angle map.")


class QIAFI(base.BaseCommand):
    """
    help for myInterface

    Example 1
    -------
    >>> from qipype.interfaces.relax import QiAfi
    >>> interface = QiAfi(prefix='nipype_', param_file='spgr_params.json')

    """

    _cmd = 'qi afi'
    input_spec = QIAFIInputSpec
    output_spec = QIAFIOutputSpec

    # If the command requires a json input file
    def _format_arg(self, name, spec, value):
        return self._process_params(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['b1_rel_map'] = os.path.abspath(
            self._add_prefix('AFI_B1.nii.gz'))
        outputs['b1_act_map'] = os.path.abspath(
            self._add_prefix('AFI_angle.nii.gz'))
        return outputs
