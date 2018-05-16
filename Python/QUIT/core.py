#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT CoreProgs.

Implemented:
    - qinewimage
    
Requires that the QUIT tools are in your your system path

By: Emil Ljungberg and Tobias Wood
"""

from __future__ import (print_function, division, unicode_literals,
                        absolute_import)

from nipype.interfaces.base import CommandLineInputSpec, CommandLine, TraitedSpec, File, traits, isdefined
import json
import os
from .base import QUITCommand, QUITCommandInputSpec

############################ qinewimage ############################

class QiNewImageInputSpec(QUITCommandInputSpec):
        
    # Options
    ndims = traits.Int(desc='Image dimension, default 3', argstr='--dims=%d')
    imsize = traits.Int(desc='Image size', argstr='--size=%d')
    voxel_spacing = traits.Float(desc='Voxel spacing', argstr='--spacing=%f')
    origin = traits.Float(desc='Image origin', argstr='--origin=%f')
    im_fill = traits.Float(desc='Fill with value', argstr='--fill=%f')
    im_gradient = traits.String(desc='Fill with value (dim, low, high)', argstr='--step=%s')
    im_steps = traits.String(desc='Fill with discrete steps (dim, low,high, steps)', argstr='--step=%s')
    wrap = traits.Float(desc='Wrap image values at the given value', argstr='--wrap=%f')

    # Output file
    out_file = traits.File(desc='Output file', exists=False, position=-1, argstr='%s', mandatory=True)

class QiNewImageOutputSpec(TraitedSpec):
    out_file = File(desc="Simulated Image")

class QiNewImage(CommandLine):
    """
    Produce a new image with qinewimage

    Example usage
    -------
    >>> from QUIT.nipype.CoreProgs import QiNewImage
    >>> qinewimage = QiNewImage(out_file='test.nii', imsize=256)
    >>> sim_res = qinewimage.run()
    >>> print(sim_res.outputs)
    """

    _cmd = 'qinewimage'
    input_spec = QiNewImageInputSpec
    output_spec = QiNewImageOutputSpec


    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = os.path.abspath('out_file.nii.gz')
        
        return outputs