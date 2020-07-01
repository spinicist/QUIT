#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT MT module
Requires that the QUIT tools are in your your system path
"""

from os import path
from nipype.interfaces.base import TraitedSpec, DynamicTraitedSpec, File, traits, isdefined
from .. import base

MTSat, MTSatSim = base.Command('MTSat', 'qi mtsat', 'MTSat',
                               varying=['PD', 'R1', 'delta'],
                               fixed=['B1'],
                               files=['PDw', 'T1w', 'MTw'])

qMT, qMTSim = base.Command('qMT', 'qi qmt', 'QMT',
                           varying=['M0_f', 'F_over_R1_f',
                                    'T2_b', 'T1_f_over_T2_f', 'k'],
                           derived=['PD', 'T1_f', 'T2_f',
                                    'T2_b', 'k_bf', 'f_b'],
                           fixed=['f0', 'B1', 'T1'],
                           extra={'lineshape': traits.String(argstr='--lineshape=%s', mandatory=True,
                                                             desc='Gauss/Lorentzian/SuperLorentzian/path to JSON file')})


eMT, eMTSim = base.Command('eMT', 'qi ssfp_emt', 'EMT',
                           varying=['PD', 'f_b', 'k_bf', 'T1_f', 'T2_f'],
                           fixed=['f0', 'B1'],
                           files=['G', 'a', 'b'],
                           extra={'T2b': traits.Float(desc='T2 of bound pool', argstr='--T2b=%f')})


def Lorentzian(name, pools):
    ###
    #  Lorentzian fitting is a bit special as it can have a variable number of pools, and the
    #  varying parameters will change name depending on the pools. Hence provide a function
    #  to create the correct Fit/Sim types for the specified pools
    #  The specified name must match the name that is used for the class in the module you call this
    #  function from if you ever intend to pickle it (which nipype does)
    ###
    """
    Returns a type for fitting a Lorentzian function to a Z-spectrum with the specified pools
    """
    pars = []
    for pool in pools:
        poolname = pool['name']
        for par in ['f0', 'fwhm', 'A']:
            pname = '{}_{}'.format(poolname, par)
            pars.append(pname)

    def T_init(self, **kwargs):
        base.FitCommand.__init__(self, **kwargs)
        self._json['pools'] = pools

    Fit, Sim = base.Command(name, 'qi lorentzian --pools={}'.format(len(pools)), 'LTZ',
                            varying=pars,
                            extra={'additive': traits.Bool(argstr='--add', desc='Use an additive instead of subtractive model'),
                                   'Zref': traits.Float(argstr='--zref=%f', desc='Set reference Z-spectrum value (usually 1 or 0)')},
                            init=T_init)
    return Fit, Sim


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

############################ qi mtr ############################


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
