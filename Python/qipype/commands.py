#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT MT module
Requires that the QUIT tools are in your your system path
"""

from os import path
from nipype.interfaces.base import TraitedSpec, DynamicTraitedSpec, File, traits, isdefined
from . import base

############################ qinewimage ############################


class NewImageInputSpec(base.InputBaseSpec):
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


class NewImage(base.CommandLine):
    """
    Produce a new image with qinewimage

    Example usage
    -------
    >>> from quit.nipype.CoreProgs import QINewImage
    >>> qinewimage = QINewImage(out_file='test.nii', imsize=256)
    >>> sim_res = qinewimage.run()
    >>> print(sim_res.outputs)
    """

    _cmd = 'qi newimage'
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


class DiffInputSpec(base.InputBaseSpec):

    # Options
    in_file = File(desc='Input file', argstr='--input=%s',
                   exists=True, mandatory=True)
    baseline = File(desc='Baseline file', argstr='--baseline=%s',
                    exists=True, mandatory=True)
    noise = traits.Float(desc='Added noise level', argstr='--noise=%f')
    abs_diff = traits.Bool(
        desc='Use absolute difference, not relative', argstr='--abs')


class DiffOutputSpec(TraitedSpec):
    out_diff = traits.Float(desc='Image difference')


class Diff(base.CommandLine):
    """
    Compare two images
    """

    _cmd = 'qi diff'
    input_spec = DiffInputSpec
    output_spec = DiffOutputSpec
    terminal_output = 'allatonce'

    def aggregate_outputs(self, runtime=None, needed_outputs=None):
        outputs = self._outputs()
        outputs.out_diff = float(runtime.stdout)
        return outputs

############################ NOISE ESTIMATION ############################


class NoiseInputSpec(base.InputBaseSpec):
    in_file = File(desc='Input file', argstr='%s',
                   position=0, mandatory=True, exists=True)
    mask_file = File(desc='Noise mark', argstr='--mask=%s', exists=True)
    region = traits.String(
        desc='Noise region. Argument should be a string "start_x,start_y,start_z,size_x,size_y,size_z"', argstr='--region=%s')
    meansqr = traits.Bool(
        desc='Return mean of squared values, not standard deviation', argstr='--meansqr')


class NoiseOutputSpec(TraitedSpec):
    noise = traits.Float(desc='Measured noise')


class Noise(base.CommandLine):
    """
    Calculate noise levels in an image
    """

    _cmd = 'qi noise_est'
    input_spec = NoiseInputSpec
    output_spec = NoiseOutputSpec
    terminal_output = 'allatonce'

    def aggregate_outputs(self, runtime=None, needed_outputs=None):
        outputs = self._outputs()
        outputs.noise = float(runtime.stdout)
        return outputs


############################ qidream ############################


class DreamInputSpec(base.InputSpec):
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


class DreamOutputSpec(TraitedSpec):
    b1_rel_map = File(desc="The relative flip-angle map.")
    b1_act_map = File(desc="The actual achieved angle in each voxel.")


class Dream(base.BaseCommand):
    _cmd = 'qi dream'
    input_spec = DreamInputSpec
    output_spec = DreamOutputSpec

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


class AFIInputSpec(base.InputSpec):
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


class AFIOutputSpec(TraitedSpec):
    # Specify which outputs there are
    b1_rel_map = File(desc="The relative flip-angle map.")
    b1_act_map = File(desc="The actual flip-angle map.")


class AFI(base.BaseCommand):
    _cmd = 'qi afi'
    input_spec = AFIInputSpec
    output_spec = AFIOutputSpec

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

############################ B1- via Papp Method ############################
# Implemented but not tested #


class B1MinusInputSpec(base.InputSpec):
    in_file = File(exists=True, argstr='%s', mandatory=True,
                   position=0, desc='Input file')


class B1MinusOutputSpec(TraitedSpec):
    out_file = File('B1minus.nii.gz',
                    desc="The B1 minus map.", usedefault=True)


class B1Minus(base.BaseCommand):
    _cmd = 'qi b1_papp'
    input_spec = B1MinusInputSpec
    output_spec = B1MinusOutputSpec


############################ FieldMap ############################


class FieldmapInputSpec(base.InputSpec):
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


class Fieldmap(base.BaseCommand):
    """
    Fieldmap via complex division

    """

    _cmd = 'qi fieldmap'
    input_spec = FieldmapInputSpec
    output_spec = FieldmapOutputSpec

    def _format_arg(self, name, spec, value):
        return self._process_params(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['fieldmap'] = os.path.abspath(
            self._add_prefix('Fieldmap.nii.gz'))

        return outputs

############################ qi_lineshape ############################


class LineshapeInputSpec(base.InputBaseSpec):
    # Options
    out_file = traits.File(argstr='%s', mandatory=True,
                           exists=False, position=-1, desc='Output file')
    lineshape = traits.String(argstr='--lineshape=%s', mandatory=True,
                              desc='Gauss/Lorentzian/SuperLorentzian')
    t2b = traits.Float(argstr='--T2b=%f',
                       desc='Nominal T2 of bound-pool (default 10µs)')
    frq_count = traits.Int(argstr='--frq_count=%d',
                           desc='Number of frequencies in table (default 10)')
    frq_start = traits.Float(argstr='--frq_start=%f',
                             desc='Start frequency for table (default 1 kHZ)')
    frq_space = traits.Float(
        argstr='--frq_space=%f', desc='Spacing of frequencies in table (default 1kHz)')


class LineshapeOutputSpec(TraitedSpec):
    out_file = File(desc='JSON Lineshape file')


class Lineshape(base.BaseCommand):
    """
    Pre-calculate lineshapes and write them out for use with qMT
    """

    _cmd = 'qi lineshape'
    input_spec = LineshapeInputSpec
    output_spec = LineshapeOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = path.abspath(self.inputs.out_file)
        return outputs

############################ qi_zspec ############################


class ZSpecInputSpec(base.InputSpec):
    # Input nifti
    in_file = File(exists=True, argstr='%s', mandatory=True,
                   position=-1, desc='Input file')
    # Options
    f0_map = File(
        desc='Fieldmap (in same units as frequencies)', argstr='--f0=%s')
    ref = File(desc='Reference image for image output', argstr='--ref=%s')
    asym = traits.Bool(
        desc='Output MT-asymmetry spectrum instead of Z-spectrum', argstr='--asym')
    order = traits.Int(
        desc='Interpolation order (default 3)', argstr='--order=%d')
    out_file = traits.String(
        desc='Output filename (default is input_interp)', argstr='--out=%s')


class ZSpecOutputSpec(TraitedSpec):
    out_file = File(desc="Path to interpolated Z-spectrum/MTA-spectrum")


class ZSpec(base.BaseCommand):
    """
    Interpolate a Z-spectrum (with correction for off-resonance)

    """

    _cmd = 'qi zspec_interp'
    input_spec = ZSpecInputSpec
    output_spec = ZSpecOutputSpec

    def __init__(self, in_freqs=[], out_freqs=[], **kwargs):
        super().__init__(**kwargs)
        self._json = {'input_freqs': in_freqs, 'output_freqs': out_freqs}

    def _list_outputs(self):
        outputs = self.output_spec().get()
        if isdefined(self.inputs.out_file):
            fname = self._gen_fname(self.inputs.out_file)
        else:
            fname = self._gen_fname(self.inputs.in_file, suffix='_interp')
        outputs['out_file'] = path.abspath(fname)
        return outputs

############################ MTR ############################


class MTRInputSpec(base.InputSpec):
    # Options
    in_file = File(exists=True, argstr='%s', mandatory=True,
                   position=-1, desc='Input file')


class MTROutputSpec(DynamicTraitedSpec):
    pass


class MTR(base.BaseCommand):
    """
    Calculate Magnetization Transfer Ratios
    """

    _cmd = 'qi mtr'
    input_spec = MTRInputSpec
    output_spec = MTROutputSpec

    def __init__(self, contrasts={}, **kwargs):
        super(MTR, self).__init__(**kwargs)
        if contrasts:
            self._json = {'contrasts': contrasts}
            for con in contrasts:
                cn = con['name']
                setattr(self.output_spec, cn, File('%s.nii.gz' %
                                                   cn, desc='Path to %s' % cn, usedefault=True))
        else:
            setattr(self.output_spec, 'MTR', File('MTR.nii.gz',
                                                  desc='Path to MTR file', usedefault=True))

    def _list_outputs(self):
        outputs = self.output_spec().get()
        if self._json:
            for con in self._json['contrasts']:
                cn = con['name']
                outputs[cn] = '%s.nii.gz' % cn
        else:
            outputs['MTR'] = 'MTR.nii.gz'
        return self._add_prefixes(outputs)

############################ Z-SHIM ############################


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

############################ MP2RAGE ############################


class MP2RAGEInputSpec(base.InputSpec):
    # Inputs
    in_file = File(exists=True, argstr='%s', mandatory=True,
                   position=0, desc='Path to complex MP-RAGE data')

    # Commonly used options
    threads = traits.Int(
        desc='Use N threads (default=hardware)', argstr='--threads=%d')
    prefix = traits.String(
        desc='Add a prefix to output filenames', argstr='--out=%s')
    beta = traits.Float(desc='Regularisation paramter', argstr='--beta=%f')


class MP2RAGEOutputSpec(TraitedSpec):
    # Specify which outputs there are
    uni_file = File('MP2_UNI.nii.gz',
                    desc='The Uniform MP2 contrast image', usedefault=True)
    T1_map = File('MP2_T1.nii.gz', desc='T1 Map', usedefault=True)


class MP2RAGE(base.BaseCommand):
    """
    Interface for qimp2rage

    Example 1
    -------
    >>> from qipype.interfaces.relax import QIMP2RAGE
    >>> interface = QIMP2RAGE(prefix='nipype_', param_file='spgr_params.json')

    """

    _cmd = 'qi mp2rage'
    input_spec = MP2RAGEInputSpec
    output_spec = MP2RAGEOutputSpec
