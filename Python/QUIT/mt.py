#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT MT module
Requires that the QUIT tools are in your your system path
"""

from os import path
from nipype.interfaces.base import TraitedSpec, DynamicTraitedSpec, File, traits, isdefined
from . import base as QI

############################ qi_lineshape ############################


class LineshapeInputSpec(QI.InputBaseSpec):
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


class Lineshape(QI.BaseCommand):
    """
    Pre-calculate lineshapes and write them out for use with qMT
    """

    _cmd = 'qi_lineshape'
    input_spec = LineshapeInputSpec
    output_spec = LineshapeOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = path.abspath(self.inputs.out_file)
        return outputs

############################ qi_lorentzian ############################


class LorentzianInputSpec(QI.FitInputSpec):
    # Options
    pools = traits.Int(argstr='--pools=%d',
                       desc='Number of pools, set from pool arg')
    additive = traits.Bool(argstr='--add',
                           desc='Use an additive instead of subtractive model')
    Zref = traits.Float(argstr='--zref=%f',
                        desc='Set reference Z-spectrum value (usually 1 or 0)')


class LorentzianOutputSpec(DynamicTraitedSpec):
    residual_map = File('LTZ_SoS_residual.nii.gz',
                        desc='Path to residual map', usedefault=True)


class Lorentzian(QI.FitCommand):
    """
    Fit a Lorentzian function to a Z-spectrum
    """

    _cmd = 'qi_lorentzian'
    input_spec = LorentzianInputSpec
    output_spec = LorentzianOutputSpec

    def __init__(self, pools={}, **kwargs):
        super(Lorentzian, self).__init__(pools=len(pools), **kwargs)
        self._json['pools'] = pools
        for pool in pools:
            pn = pool['name']
            setattr(self.output_spec, pn + '_f0', File('LTZ_%s_f0.nii.gz' %
                                                       pn, desc='Path to %s Δf0 map' % pn, usedefault=True))
            setattr(self.output_spec, pn + '_fwhm', File('LTZ_%s_fwhm.nii.gz' %
                                                         pn, desc='Path to %s FWHM map' % pn, usedefault=True))
            setattr(self.output_spec, pn + '_A', File('LTZ_%s_A.nii.gz' %
                                                      pn, desc='Path to %s A map' % pn, usedefault=True))

    def _list_outputs(self):
        outputs = self.output_spec().get()
        for pool in self._json['pools']:
            pn = pool['name']
            outputs[pn + '_f0'] = 'LTZ_%s_f0.nii.gz' % pn
            outputs[pn + '_fwhm'] = 'LTZ_%s_fwhm.nii.gz' % pn
            outputs[pn + '_A'] = 'LTZ_%s_A.nii.gz' % pn
        return self._add_prefixes(outputs)


class LorentzianSimInputSpec(QI.SimInputSpec):
    # Options
    # Options
    pools = traits.Int(argstr='--pools=%d',
                       desc='Number of pools, set from pool arg')
    additive = traits.Bool(argstr='--add',
                           desc='Use an additive instead of subtractive model')
    Zref = traits.Float(argstr='--zref=%f',
                        desc='Set reference Z-spectrum value (usually 1 or 0)')


class LorentzianSim(QI.SimCommand):
    _cmd = 'qi_lorentzian'
    input_spec = LorentzianSimInputSpec
    output_spec = QI.SimOutputSpec

    def __init__(self, pools={}, **kwargs):
        self._param_files = []
        for pool in pools:
            pn = pool['name']
            self._param_files.append(pn + '_f0')
            self._param_files.append(pn + '_fwhm')
            self._param_files.append(pn + '_A')
        super().__init__(pools=len(pools), **kwargs)
        self._json['pools'] = pools

############################ qi_qmt ############################


class qMTInputSpec(QI.FitInputSpec):
    # Inputs
    t1_map = File(exists=True, argstr='%s', mandatory=True,
                  position=-2, desc='Path to T1 map')

    # Options
    lineshape = traits.String(argstr='--lineshape=%s', mandatory=True,
                              desc='Gauss/Lorentzian/SuperLorentzian/path to JSON file')
    f0_map = File(desc='f0 map (Hertz)', argstr='--f0=%s', exists=True)
    b1_map = File(desc='B1 map (ratio) file', argstr='--B1=%s')


class qMTOutputSpec(TraitedSpec):
    pd_map = File('QMT_PD.nii.gz', desc="Path to PD map", usedefault=True)
    T1_f_map = File('QMT_T1_f.nii.gz',
                    desc="Path to T1 of free pool", usedefault=True)
    T2_f_map = File('QMT_T2_f.nii.gz',
                    desc="Path to T2 of free pool", usedefault=True)
    T2_b_map = File('QMT_T2_b.nii.gz',
                    desc="Path to T2 of bound pool", usedefault=True)
    k_bf_map = File('QMT_k_bf.nii.gz',
                    desc="Path to exchange rate from bound to free pool", usedefault=True)
    f_b_map = File('QMT_f_b.nii.gz',
                   desc="Path to bound pool fraction", usedefault=True)
    residual_map = File('QMT_SoS_residual.nii.gz',
                        desc="Path to residual map", usedefault=True)


class qMT(QI.FitCommand):
    """
    Fit the Ramani model to a Z-spectrum
    """

    _cmd = 'qi_qmt'
    input_spec = qMTInputSpec
    output_spec = qMTOutputSpec


class qMTSimInputSpec(QI.SimInputSpec):
    # Inputs
    t1_map = File(argstr='%s', mandatory=True,
                  position=-2, desc='Path to T1 map')

    # Options
    lineshape = traits.String(argstr='--lineshape=%s', mandatory=True,
                              desc='Gauss/Lorentzian/SuperLorentzian/path to JSON file')
    f0_map = File(desc='f0 map (Hertz)', argstr='--f0=%s', exists=True)
    b1_map = File(desc='B1 map (ratio) file', argstr='--B1=%s')


class qMTSim(QI.SimCommand):
    _cmd = 'qi_qmt'
    _param_files = ['PD', 'T1_f', 'T2_f', 'T2_b', 'k_bf', 'f_b']
    input_spec = qMTSimInputSpec
    output_spec = QI.SimOutputSpec

############################ qi_zspec ############################


class ZSpecInputSpec(QI.InputSpec):
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


class ZSpec(QI.BaseCommand):
    """
    Interpolate a Z-spectrum (with correction for off-resonance)

    """

    _cmd = 'qi_zspec_interp'
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
