#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT CoreProgs.

Implemented:
    - qinewimage
    
Requires that the QUIT tools are in your your system path

By: Emil Ljungberg and Tobias Wood
"""

from nipype.interfaces.base import CommandLine, TraitedSpec, File, traits
from os import path
from . import base as QI

############################ qinewimage ############################


class NewImageInputSpec(QI.InputBaseSpec):
    # Options
    img_size = traits.List(minsize=2, maxsize=4, mandatory=True,
                           desc='Image size', argstr='--size=%s', sep=',')
    voxel_spacing = traits.Float(desc='Voxel spacing', argstr='--spacing=%f')
    origin = traits.Float(desc='Image origin', argstr='--origin=%f')
    fill = traits.Float(desc='Fill with value', argstr='--fill=%f')
    grad_dim = traits.Int(
        desc='Fill with gradient along dimension', argstr='--grad_dim=%d')
    grad_vals = traits.Tuple(
        desc='Gradient start/end values', argstr='--grad_vals=%f,%f')
    grad_steps = traits.Int(
        desc='Gradient in N discrete steps', argstr='--steps=%s')
    wrap = traits.Float(
        desc='Wrap image values at the given value', argstr='--wrap=%f')

    # Output file
    out_file = traits.File(desc='Output file', exists=False,
                           position=-1, argstr='%s', mandatory=True)


class NewImageOutputSpec(TraitedSpec):
    out_file = File(desc="Simulated Image")


class NewImage(CommandLine):
    """
    Produce a new image with qinewimage

    Example usage
    -------
    >>> from QUIT.nipype.CoreProgs import QINewImage
    >>> qinewimage = QINewImage(out_file='test.nii', imsize=256)
    >>> sim_res = qinewimage.run()
    >>> print(sim_res.outputs)
    """

    _cmd = 'qinewimage'
    input_spec = NewImageInputSpec
    output_spec = NewImageOutputSpec

    def _parse_inputs(self, skip=None):
        dim_arg = '--dims=%d' % len(self.inputs.img_size)
        return [dim_arg, ] + super()._parse_inputs(skip)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = path.abspath(self.inputs.out_file)
        return outputs

############################ qidiff ############################


class DiffInputSpec(QI.InputBaseSpec):

    # Options
    in_file = traits.String(desc='Input file', argstr='--input=%s')
    baseline = traits.String(desc='Baseline file', argstr='--baseline=%s')
    noise = traits.Float(desc='Added noise level', argstr='--noise=%f')
    abs_diff = traits.Bool(
        desc='Use absolute difference, not relative', argstr='--abs')


class DiffOutputSpec(TraitedSpec):
    out_diff = traits.Float(desc='Image difference')


class Diff(CommandLine):
    """
    Compare two images
    """

    _cmd = 'qidiff'
    input_spec = DiffInputSpec
    output_spec = DiffOutputSpec

    def aggregate_outputs(self, runtime=None, needed_outputs=None):
        outputs = self._outputs()
        outputs.out_diff = float(runtime.stdout)
        return outputs
