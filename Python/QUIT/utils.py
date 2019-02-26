#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT utilities.

Implemented:
    - qiaffine
    - qimask
    - qi_rfprofile
    
To be implemented:
    - qi_coil_combine
    - qicomplex
    - qihdr
    - qikfilter
    - qipolyfit
    - qipolyimg

Requires that the QUIT tools are in your your system path
"""

from __future__ import (print_function, division, unicode_literals,
                        absolute_import)

from nipype.interfaces.base import CommandLineInputSpec, CommandLine, TraitedSpec, File, traits, isdefined
# from quit_nipype_utils import parse_param_dict
from .base import QUITCommand, QUITCommandInputSpec

import json
import os

############################ qi_rfprofile ############################


class QiRFprofileInputSpec(QUITCommandInputSpec):

    # Inputs
    in_file = File(exists=True, argstr='%s', mandatory=True,
                   position=0, desc='Input B1+ file')
    param_file = File(desc='Parameter .json file', position=2, argstr='--json=%s',
                      xor=['param_dict'], mandatory=True, exists=True)
    param_dict = traits.Dict(desc='Parameter dictionary', position=2,
                             argstr='', mandatory=True, xor=['param_file'])

    # Outputs
    out_file = File(exists=False, argstr='%s', mandatory=True,
                    desc='Input B1+ file', position=1)

    # Options
    mask_file = File(
        desc='Only process voxels within the mask', argstr='--mask=%s')

    # Commonly used options
    debug = traits.Bool(desc='Output debugging messages', argstr='-v')
    threads = traits.Int(
        desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(
        desc='Add a prefix to output filenames', argstr='--out=%s')


class QiRFprofileOutputSpec(TraitedSpec):
    # Specify which outputs there are
    out_file = File(desc="output relative B1 file")


class QiRFprofile(QUITCommand):
    """
    help for myInterface

    Example 1
    -------
    >>> from QUIT.nipype.utils import QiRFprofile
    >>> interface = QiRFprofile()

    """

    _cmd = 'qi_rfprofile'
    input_spec = QiRFprofileInputSpec
    output_spec = QiRFprofileOutputSpec

    def _format_arg(self, name, spec, value):
        return self._process_params(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()

        outputs['out_file'] = os.path.abspath(
            self._add_prefix(self.inputs.out_file))

        return outputs

############################ qiaffine ############################


class QiAffineInputSpec(QUITCommandInputSpec):

    # Inputs
    in_file = File(exists=True, argstr='%s', mandatory=True,
                   position=0, desc='Source File')

    # Outputs
    out_file = File(exists=False, argstr='%s', mandatory=True,
                    position=1, desc='Destination File')

    # Options
    xfm_file = File(exists=False, argstr='--tfm=%s', mandatory=True,
                    desc='Write out the transformation to a file')
    offX = traits.Int(desc='Translate origin in X direction',
                      argstr='-offX=%d')
    offY = traits.Int(desc='Translate origin in Y direction',
                      argstr='-offY=%d')
    offZ = traits.Int(desc='Translate origin in Z direction',
                      argstr='-offZ=%d')
    rotX = traits.Int(
        desc='Rotate about X-axis by angle (degrees)', argstr='-rotX=%d')
    rotY = traits.Int(
        desc='Rotate about Y-axis by angle (degrees)', argstr='-rotY=%d')
    rotZ = traits.Int(
        desc='Rotate about Z-axis by angle (degrees)', argstr='-rotZ=%d')
    scale = traits.Float(desc='Scale by a constant', argstr='--scale=%f')
    permute = traits.String(
        desc='Permute axes, e.g. 2,0,1. Negative values mean flip as well', argstr='--permute=%s')
    flip = traits.String(
        desc='Flip an axis, e.g. 0,1,0. Occurs AFTER any permutation.', argstr='--flip=%s')

    # This should probably be enum instead
    center = traits.String(
        desc='Set the origin to geometric center (geo) or (cog)', argstr='--center=%s')


class QiAffineOutputSpec(TraitedSpec):
    # Specify which outputs there are
    out_file = File(desc="output relative B1 file")


class QiAffine(QUITCommand):
    """
    help for myInterface

    Example 1
    -------
    >>> from QUIT.nipype.utils import QiRFprofile
    >>> interface = QiRFprofile()

    """

    _cmd = 'qi_rfprofile'
    input_spec = QiRFprofileInputSpec
    output_spec = QiRFprofileOutputSpec

    def _format_arg(self, name, spec, value):
        return self._process_params(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = os.path.abspath(
            self._add_prefix(self.inputs.out_file))
        return outputs

############################ qimask ############################


class QiMaskInputSpec(QUITCommandInputSpec):
    # Inputs
    in_file = File(exists=True, argstr='%s', mandatory=True,
                   position=0, desc='Input File')

    # Outputs
    out_file = File(exists=False, argstr='%s', mandatory=True,
                    position=1, desc='Set output filename, default is input + _mask')

    # Options
    vol = traits.Int(desc='Choose volume to mask in multi-volume file. Default 1, -1 selects last volume',
                     argstr='--volume=%d')
    complex_data = traits.Bool(
        desc='Input data is complex, take magnitude first', argstr='--complex')
    lower = traits.Float(desc="Specify lower intensity threshold for 1st stage, otherwise Otsu's method is used",
                         argstr='--lower=%f')
    upper = traits.Float(desc="Specify upper intensity threshold for 1st stage, otherwise Otsu's method is used",
                         argstr='--upper=%f')
    rats = traits.Float(desc="Perform the RATS step, argument is size threshold for connected component",
                        argstr='--rats=%f')
    fill_holes = traits.Int(
        desc="Fill holes in thresholded mask with radius N", argstr='--fillh=%d')


class QiMaskOutputSpec(TraitedSpec):
    # Specify which outputs there are
    mask_file = File(desc="Output mask file")


class QiMask(CommandLine):
    """
    help for myInterface

    Example 1
    -------
    >>> from QUIT.nipype.utils import QiRFprofile
    >>> interface = QiRFprofile()

    """

    _cmd = 'qimask'
    input_spec = QiMaskInputSpec
    output_spec = QiMaskOutputSpec

    def _format_arg(self, name, spec, value):
        return self._process_params(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()

        # Specify output files
        if self.inputs.out_file:
            output_fname = self.inputs.out_file
        else:
            input_bname = os.path.basename(self.inputs.in_file)
            input_fname = os.path.splitext(os.path.split(input_bname)[-1])
            output_fname = input_fname + '_mask.nii.gz'

        outputs['mask_file'] = os.path.abspath(output_fname)

        return outputs


############################ qi_coil_combine ############################
# < To be implemented > #

############################ qicomplex ############################
class ComplexInputSpec(QUITCommandInputSpec):

    # Options
    mag = traits.String(desc='Magnitude input', argstr='--mag=%s', exists=True)
    pha = traits.String(desc='Phase input', argstr='--pha=%s', exists=True)
    real = traits.String(desc='Real input', argstr='--real=%s', exists=True)
    imag = traits.String(desc='Imaginary input',
                         argstr='--imag=%s', exists=True)
    x = traits.String(desc='Complex input', argstr='--complex=%s', exists=True)
    realimag = traits.String(
        desc='Real/Imaginary input', argstr='--realimag=%s', mandatory=True, exists=True)

    complex_out_file = File(argstr='--COMPLEX=%s',
                            desc='Output complex file',
                            genfile=True,
                            hash_files=False)

    fixge = traits.Bool(
        desc='Fix GE FFT-shift bug (negate alternate slices)', argstr='--fixge')
    negate = traits.Bool(desc='Multiply by -1', argstr='--negate')
    conjugate = traits.Bool(desc='Conjugate data', argstr='--conjugate')


class ComplexOutputSpec(TraitedSpec):
    # Specify which outputs there are
    complex_out_file = File()


class Complex(QUITCommand):
    """
    Deals with magnitude/phase/real/imaginary/complex data

    Example 1
    -------
    >>> from QUIT.nipype.utils import Complex
    >>> interface = Complex()

    """

    _cmd = 'qicomplex'
    input_spec = ComplexInputSpec
    output_spec = ComplexOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        inputs = self.input_spec().get()

        if not isdefined(inputs['complex_out_file']):
            inputs['complex_out_file'] = self._gen_fname(
                self.inputs.realimag, suffix='_x')
            outputs['complex_out_file'] = inputs['complex_out_file']
        return outputs

    def _gen_filename(self, name):
        if name == 'complex_out_file':
            return self._list_outputs()[name]
        return None

############################ qihdr ############################
# < To be implemented > #

############################ qikfilter ############################
# < To be implemented > #

############################ qipolyfit ############################
# < To be implemented > #

############################ qipolyimg ############################
# < To be implemented > #
