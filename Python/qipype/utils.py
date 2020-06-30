#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT utilities.

Requires that the QUIT tools are in your your system path
"""

from json import dump, loads
from os import path
from nipype.interfaces.base import CommandLine, TraitedSpec, File, traits, isdefined
from . import base

############################### qi_pca ###############################


class PCAInputSpec(base.InputSpec):
    in_file = File(argstr='%s', mandatory=True, exists=True,
                   position=-1, desc='Input file for PCA denoising')
    retain = traits.Int(argstr='--retain=%d', desc='Number of PCs to retain')
    mask_file = File(argstr='--mask=%s', exists=True, desc='Mask file')
    projections_file = traits.String(
        argstr='--project=%s', desc='File to save projections to')
    pc_json_file = traits.String(
        argstr='--save_pcs=%s', desc='Save Principal Components to JSON file')
    out_file = traits.String(
        argstr='--out=%s', desc='Name of output file (default is input_pca)')


class PCAOutputSpec(TraitedSpec):
    out_file = File(desc='Denoised image')
    projections_file = File(desc='Projections file')
    pc_json_file = File(desc='JSON containing Principal Components')


class PCA(base.BaseCommand):
    """
    Denoise an image using Principal Components
    """
    _cmd = 'qi pca'
    input_spec = PCAInputSpec
    output_spec = PCAOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        if isdefined(self.inputs.out_file):
            outputs['out_file'] = path.abspath(self.inputs.out_file)
        else:
            p, f = path.split(self.inputs.in_file)
            fname, ext = path.splitext(f)
            if ext == '.gz':
                fname = path.splitext(fname)[0]
            outputs['out_file'] = path.abspath(fname + '_pca.nii.gz')
        if isdefined(self.inputs.pc_json_file):
            outputs['pc_json_file'] = path.abspath(self.inputs.pc_json_file)
        if isdefined(self.inputs.projections_file):
            outputs['projections_file'] = path.abspath(
                self.inputs.projections_file)
        return outputs


############################ qi_rfprofile ############################


class RFProfileInputSpec(base.InputSpec):
    rf = traits.Dict(desc='Dictionary with rf_pos and rf_vals lists', argstr='',
                     mandatory=True)

    # Input nifti
    in_file = File(exists=True, argstr='%s', mandatory=True,
                   position=-2, desc='Input B1+ map')
    out_file = File(argstr='%s', mandatory=True,
                    position=-1, desc='Output slab profile')

    threads = traits.Int(
        desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(
        desc='Add a prefix to output filenames', argstr='--out=%s')
    mask_file = File(
        desc='Only process voxels within the mask', argstr='--mask=%s')
    center = traits.List(
        desc='Set center point of slab profile to center of mask', argsstr='--center')


class RFProfileOutputSpec(TraitedSpec):
    # Specify which outputs there are
    out_file = File(desc="output slab profile file")


class RFProfile(base.BaseCommand):
    """
    help for myInterface

    Example 1
    -------
    >>> from quit.nipype.utils import RFProfile
    >>> interface = RFProfile()

    """

    _cmd = 'qi rfprofile'
    input_spec = RFProfileInputSpec
    output_spec = RFProfileOutputSpec

    def _parse_inputs(self, skip=None):
        # Make sequence dictionary into a .json file for input to interface
        if isdefined(self.inputs.rf):
            self._json = self.inputs.rf
        return super()._parse_inputs(skip)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = path.abspath(self.inputs.out_file)
        return outputs

############################ qiaffine ############################


class AffineInputSpec(base.InputSpec):

    # Inputs
    in_file = File(exists=True, argstr='%s', mandatory=True,
                   position=0, desc='Source File')

    # Outputs
    out_file = File(exists=False, argstr='%s', mandatory=True,
                    position=1, desc='Destination File')

    # Options
    xfm_file = File(exists=False, argstr='--tfm=%s',
                    desc='Write out the transformation to a file')
    translate = traits.String(desc='Translations in "X,Y,Z" directions',
                              argstr='--trans=%s')
    rotate = traits.String(
        desc='Rotations about "X,Y,Z" in degrees', argstr='--rotate=%s')
    scale = traits.Float(desc='Scale by a constant', argstr='--scale=%f')
    permute = traits.String(
        desc='Permute axes, e.g. 2,0,1. Negative values mean flip as well', argstr='--permute=%s')
    flip = traits.String(
        desc='Flip an axis, e.g. 0,1,0. Occurs AFTER any permutation.', argstr='--flip=%s')

    # This should probably be enum instead
    center = traits.String(
        desc='Set the origin to geometric center (geo) or (cog)', argstr='--center=%s')


class AffineOutputSpec(TraitedSpec):
    # Specify which outputs there are
    out_file = File(desc="Transformed file")


class Affine(CommandLine):
    """
    help for myInterface

    Example 1
    -------
    >>> from quit.nipype.utils import Affine
    >>> interface = Affine()

    """

    _cmd = 'qi affine'
    input_spec = AffineInputSpec
    output_spec = AffineOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = path.abspath(self.inputs.out_file)
        return outputs

############################ qimask ############################


class MaskInputSpec(base.InputBaseSpec):
    # Inputs
    in_file = File(exists=True, argstr='%s', mandatory=True,
                   position=0, desc='Input File')

    # Outputs
    out_file = File(exists=False, argstr='--out=%s',
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
    rats_radius = traits.Int(
        desc="Starting radius for RATS (saves times), default 1", argstr='--rats_radius=%d')
    fill_holes = traits.Int(
        desc="Fill holes in thresholded mask with radius N", argstr='--fillh=%d')


class MaskOutputSpec(TraitedSpec):
    # Specify which outputs there are
    out_file = File(desc="Output mask file")


class Mask(base.BaseCommand):
    """
    help for myInterface

    Example 1
    -------
    >>> from quit.nipype.utils import RFProfile
    >>> interface = RFProfile()

    """

    _cmd = 'qi mask'
    input_spec = MaskInputSpec
    output_spec = MaskOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        if isdefined(self.inputs.out_file):
            fname = self._gen_fname(self.inputs.out_file)
        else:
            fname = self._gen_fname(self.inputs.in_file, suffix='_mask')
        outputs['out_file'] = fname
        return outputs


############################ qi_coil_combine ############################
class CoilCombineInputSpec(base.InputSpec):
    in_file = traits.File(
        desc='Input complex-valued file to coil-combine', argstr='%s', exists=True)
    composer_file = traits.File(
        desc='Short Echo-Time reference file for COMPOSER', argstr='--composer=%s', exists=True)
    hammond_coils = traits.Int(
        desc='Number of coils for Hammond method', argstr='--coils=%d')
    hammond_volume = traits.Int(
        desc='Volume to use for Hammond method', argstr='--vol=%d')
    hammond_region = traits.Str(
        desc='Region to use for Hammond method', argstr='--region=%s')


class CoilCombineOutputSpec(TraitedSpec):
    out_file = traits.File(desc='Output combined file')


class CoilCombine(base.BaseCommand):
    """
    Combine multi-coil data
    """
    _cmd = 'qi coil_combine'
    input_spec = CoilCombineInputSpec
    output_spec = CoilCombineOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        if isdefined(self.inputs.prefix):
            outputs['out_file'] = path.abspath(
                self.inputs.prefix + '_combined.nii.gz')
        else:
            p, f = path.split(self.inputs.in_file)
            fname, ext = path.splitext(f)
            if ext == '.gz':
                fname = path.splitext(fname)[0]
            outputs['out_file'] = path.abspath(fname + '_combined.nii.gz')
        return outputs

############################ qicomplex ############################


class ComplexInputSpec(base.InputSpec):

    # Options
    mag = traits.String(desc='Magnitude input', argstr='--mag=%s', exists=True)
    pha = traits.String(desc='Phase input', argstr='--pha=%s', exists=True)
    real = traits.String(desc='Real input', argstr='--real=%s', exists=True)
    imag = traits.String(desc='Imaginary input',
                         argstr='--imag=%s', exists=True)
    x = traits.String(desc='Complex input', argstr='--complex=%s', exists=True)
    realimag = traits.String(
        desc='Real/Imaginary input', argstr='--realimag=%s', exists=True)

    magnitude_out_file = File(
        argstr='--MAG=%s', desc='Output magnitude file', hash_files=False)
    real_out_file = File(
        argstr='--REAL=%s', desc='Output real file', hash_files=False)
    imag_out_file = File(
        argstr='--IMAG=%s', desc='Output imaginary file', hash_files=False)
    complex_out_file = File(argstr='--COMPLEX=%s',
                            desc='Output complex file',
                            hash_files=False)

    fix_ge = traits.Bool(
        desc='Fix GE FFT-shift bug (negate alternate slices)', argstr='--fixge')
    negate = traits.Bool(desc='Multiply by -1', argstr='--negate')
    conjugate = traits.Bool(desc='Conjugate data', argstr='--conjugate')


class ComplexOutputSpec(TraitedSpec):
    # Specify which outputs there are
    complex_out_file = File()
    magnitude_out_file = File()
    real_out_file = File()
    imag_out_file = File()


class Complex(base.BaseCommand):
    """
    Deals with magnitude/phase/real/imaginary/complex data

    Example 1
    -------
    >>> from quit.nipype.utils import Complex
    >>> interface = Complex()

    """

    _cmd = 'qi complex'
    input_spec = ComplexInputSpec
    output_spec = ComplexOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        inputs = self.inputs

        if isdefined(inputs.complex_out_file):
            outputs['complex_out_file'] = path.abspath(inputs.complex_out_file)
        if isdefined(inputs.magnitude_out_file):
            outputs['magnitude_out_file'] = path.abspath(
                inputs.magnitude_out_file)
        if isdefined(inputs.real_out_file):
            outputs['real_out_file'] = path.abspath(inputs.real_out_file)
        if isdefined(inputs.imag_out_file):
            outputs['imag_out_file'] = path.abspath(inputs.imag_out_file)
        return outputs

############################ qihdr ############################
# < To be implemented > #

############################ qikfilter ############################


class FilterInputSpec(base.InputSpec):
    in_file = File(argstr='%s', mandatory=True, exists=True,
                   position=-1, desc='Input file to fit polynomial to')
    filter_spec = traits.String(argstr='--filter=%s', mandatory=True,
                                desc='Filter to apply', multiple=True)
    filter_per_volume = traits.Bool(
        argstr='--filter_per_volume', desc='Filters are one per volume')
    complex_in = traits.Bool(argstr='--complex_in', desc='Read complex data')
    complex_out = traits.Bool(argstr='--complex_out',
                              desc='Write complex data')
    save_kernel = traits.Bool(argstr='--save_kernel',
                              desc='Save k-Space kernel')
    highpass = traits.Bool(argstr='--highpass',
                           desc='Highpass instead of lowpass')
    prefix = traits.String(
        argstr='--out=%s', desc='Output prefix (default is input filename)')
    zeropad = traits.Int(argstr='--zero_pad=%d',
                         desc='Zero-pad volume by N voxels in each direction')


class FilterOutputSpec(TraitedSpec):
    out_file = File(desc="Simulated Image")


class Filter(base.BaseCommand):
    """
    Filter an image in k-space
    """
    _cmd = 'qi kfilter'
    input_spec = FilterInputSpec
    output_spec = FilterOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        if isdefined(self.inputs.prefix):
            outputs['out_file'] = path.abspath(
                self.inputs.prefix + '_filtered.nii.gz')
        else:
            p, f = path.split(self.inputs.in_file)
            fname, ext = path.splitext(f)
            if ext == '.gz':
                fname = path.splitext(fname)[0]
            outputs['out_file'] = path.abspath(fname + '_filtered.nii.gz')
        return outputs

############################ qipolyimg ############################


class PolyImageInputSpec(base.InputSpec):
    # Options
    ref_file = File(argstr='%s', mandatory=True, exists=True, position=-2,
                    desc='Reference file for co-ordinate space and size')
    out_file = traits.File(argstr='%s', mandatory=True,
                           exists=False, position=-1, desc='Output file')
    order = traits.Int(argstr='--order=%d', mandatory=True,
                       desc='Polynomial Order')
    poly = traits.Dict(
        mandatory=True, desc='Polynomial paramters (center, scale, coeffs)', argstr='')


class PolyImageOutputSpec(TraitedSpec):
    out_file = File(desc='Output polynomial image')


class PolyImage(base.BaseCommand):
    """
    Produce a new image with qipolyimage
    """

    _cmd = 'qi polyimg'
    input_spec = PolyImageInputSpec
    output_spec = PolyImageOutputSpec

    def _format_arg(self, name, spec, value):
        """
        Make parameter dictionary into a .json file for input to interface
        """
        if name == 'poly':
            fname = '_tmp_poly.json'
            with open(fname, 'w') as outfile:
                dump(value, outfile)
            newarg = "--json=" + fname
            return newarg
        else:
            return super()._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = path.abspath(self.inputs.out_file)
        return outputs

############################ qipolyfit ############################


class PolyFitInputSpec(base.InputBaseSpec):
    in_file = File(argstr='%s', mandatory=True, exists=True,
                   position=-1, desc='Input file to fit polynomial to')
    mask_file = File(argstr='--mask=%s', exists=True,
                     desc='Mask file for fitting polynomial')
    order = traits.Int(argstr='--order=%d', mandatory=True,
                       desc='Polynomial Order')
    robust = traits.Bool(
        argstr='--robust', desc='Use robust (iterative) polynomial fit')


class PolyFitOutputSpec(TraitedSpec):
    poly = traits.Dict(desc="Polynomial parameters")


class PolyFit(base.BaseCommand):
    """
    Fit a polynomial to an image
    """
    _cmd = 'qi polyfit'
    input_spec = PolyFitInputSpec
    output_spec = PolyFitOutputSpec

    def aggregate_outputs(self, runtime=None, needed_outputs=None):
        outputs = self._outputs()
        outputs.poly = loads(runtime.stdout)
        return outputs

############################ qi_select ############################


class SelectInputSpec(base.InputBaseSpec):
    # Inputs
    in_file = File(exists=True, argstr='%s', mandatory=True,
                   position=0, desc='Input File')
    out_file = File(exists=False, argstr='%s', mandatory=True,
                    position=1, desc='Output File')
    volumes = traits.List(mandatory=True, argstr='%s', sep=',')


class SelectOutputSpec(TraitedSpec):
    # Specify which outputs there are
    out_file = File(desc="Output file")


class Select(base.BaseCommand):
    """
    Selects a set of volumes from a timeseries and writes them to another series

    """

    _cmd = 'qi select'
    input_spec = SelectInputSpec
    output_spec = SelectOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        fname = self._gen_fname(self.inputs.out_file)
        outputs['out_file'] = fname
        return outputs
