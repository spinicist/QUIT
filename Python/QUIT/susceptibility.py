#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT Susceptibility.
    
Requires that the QUIT tools are in your your system path

By: Emil Ljungberg and Tobias Wood
"""

from __future__ import (print_function, division, unicode_literals,
                        absolute_import)

from nipype.interfaces.base import CommandLineInputSpec, CommandLine, TraitedSpec, File, traits, isdefined
import json
import os
from .base import QUITCommand, QUITCommandInputSpec

############################ qi_unwrap_laplace ############################
# < To be implemented > #

############################ qi_unwrap_path ############################
# < To be implemented > #

############################ qidespot1 ############################


class FieldmapInputSpec(QUITCommandInputSpec):
    # Inputs
    input_file = File(exists=True,
                      argstr='%s',
                      mandatory=True,
                      desc='Path to input data')
    delta_te = traits.Float(
        desc='ΔTE', argstr='--delta_te=%f', mandatory=True)
    B0 = traits.Float(desc='Field strength (Tesla)', argstr='--B0=%f')

    # Options
    threads = traits.Int(
        desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(
        desc='Add a prefix to output filename', argstr='--out=%s')


class FieldmapOutputSpec(TraitedSpec):
    fieldmap = File(desc="Path to fieldmap")


class Fieldmap(QUITCommand):
    """
    Fieldmap via complex division

    """

    _cmd = 'qi_fieldmap'
    input_spec = FieldmapInputSpec
    output_spec = FieldmapOutputSpec

    def _format_arg(self, name, spec, value):
        return self._process_params(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['fieldmap'] = os.path.abspath(
            self._add_prefix('Fieldmap.nii.gz'))

        return outputs
