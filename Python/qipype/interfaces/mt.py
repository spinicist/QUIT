#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT MT module
Requires that the QUIT tools are in your your system path
"""

from os import path
from nipype.interfaces.base import TraitedSpec, DynamicTraitedSpec, File, traits, isdefined
from .. import base

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

############################ qi_lorentzian ############################


class LorentzianInputSpec(base.FitInputSpec):
    # Options
    additive = traits.Bool(argstr='--add',
                           desc='Use an additive instead of subtractive model')
    Zref = traits.Float(argstr='--zref=%f',
                        desc='Set reference Z-spectrum value (usually 1 or 0)')


def Lorentzian(pools):
    """
    Returns a type for fitting a Lorentzian function to a Z-spectrum with the specified pools
    """
    pars = []
    for pool in pools:
        poolname = pool['name']
        for par in ['f0', 'fwhm', 'A']:
            pname = '{}_{}'.format(poolname, par)
            pars.append(pname)
    ospec = base.FitOutputSpec('LTZ', pars)

    def T_init(self, **kwargs):
        base.FitCommand.__init__(self, **kwargs)
        self._json['pools'] = pools

    attrs = {'_cmd': 'qi lorentzian --pools={}'.format(len(pools)),
             'input_spec': LorentzianInputSpec,
             'output_spec': ospec,
             '__init__': T_init}

    name = 'Lorentzian{}Sim'.format(len(pools))

    T = type(name, (base.FitCommand,), attrs)
    return T


def LorentzianSim(pools):
    """
    Returns a type for simulating a Lorentzian Z-spectrum with the specified pools
    """
    pars = []
    for pool in pools:
        poolname = pool['name']
        for par in ['f0', 'fwhm', 'A']:
            pname = '{}_{}'.format(poolname, par)
            pars.append(pname)
    ispec = base.SimInputSpec('LTZ', pars, extras={'additive': traits.Bool(argstr='--add',
                                                                           desc='Use an additive instead of subtractive model'),
                                                   'Zref': traits.Float(argstr='--zref=%f',
                                                                        desc='Set reference Z-spectrum value (usually 1 or 0)')})
    ospec = base.SimOutputSpec('LTZ')

    def T_init(self, **kwargs):
        base.SimCommand.__init__(self, **kwargs)
        self._json['pools'] = pools

    attrs = {'_cmd': 'qi lorentzian --pools={}'.format(len(pools)),
             'input_spec': ispec,
             'output_spec': ospec,
             '__init__': T_init}

    name = 'Lorentzian{}Sim'.format(len(pools))

    T = type(name, (base.SimCommand,), attrs)
    return T


############################ qi_qmt ############################


class qMTInputSpec(base.FitInputSpec):
    T1_map = File(exists=True, argstr='--T1=%s',
                  mandatory=True, desc='Path to T1 map')
    f0_map = File(desc='f0 map (Hertz)', argstr='--f0=%s', exists=True)
    B1_map = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    lineshape = traits.String(argstr='--lineshape=%s', mandatory=True,
                              desc='Gauss/Lorentzian/SuperLorentzian/path to JSON file')


class qMTOutputSpec(TraitedSpec):
    pd_map = File('QMT_M0_f.nii.gz', desc="Path to M0 map", usedefault=True)
    T1_f_map = File('QMT_T1_f.nii.gz',
                    desc="Path to T1 of free pool", usedefault=True)
    T2_f_map = File('QMT_T2_f.nii.gz',
                    desc="Path to T2 of free pool", usedefault=True)
    T2_b_map = File('QMT_T2_b.nii.gz',
                    desc="Path to T2 of bound pool", usedefault=True)
    k_bf_map = File('QMT_k.nii.gz',
                    desc="Path to exchange rate", usedefault=True)
    f_b_map = File('QMT_f_b.nii.gz',
                   desc="Path to bound pool fraction", usedefault=True)
    rmse_map = File('QMT_rmse.nii.gz',
                    desc="Path to residual map", usedefault=True)


class qMT(base.FitCommand):
    """
    Fit the Ramani model to a Z-spectrum. The output parameters here are the interesting (derived)
    parameters, not the intrinsic model parameters.
    """

    _cmd = 'qi qmt'
    input_spec = qMTInputSpec
    output_spec = base.FitOutputSpec(
        'QMT', ['PD', 'T1_f', 'T2_f', 'k_bf', 'f_b'])


class qMTSim(base.SimCommand):
    """
    Simulate a Z-spectrum using the Ramani model. Note that the parameters here are the intrinsic
    model parameters, which are different to the derived parameters output from the fitting
    """
    _cmd = 'qi qmt'
    input_spec = base.SimInputSpec('qMT',
                                   varying=['M0_f', 'F_over_R1_f',
                                            'T2_b', 'T1_f_over_T2_f', 'k'],
                                   fixed=['f0', 'B1', 'T1'],
                                   extras={'lineshape': traits.String(argstr='--lineshape=%s', mandatory=True,
                                                                      desc='Gauss/Lorentzian/SuperLorentzian/path to JSON file')})
    output_spec = base.SimOutputSpec('qMT')

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

############################ qi_ssfp_emt ############################


class eMTInputSpec(base.InputSpec):
    # Inputs
    G_map = File(exists=True, argstr='%s', mandatory=True,
                 position=-3, desc='Path to G parameter map')
    a_map = File(exists=True, argstr='%s', mandatory=True,
                 position=-2, desc='Path to a parameter map')
    b_map = File(exists=True, argstr='%s', mandatory=True,
                 position=-1, desc='Path to b parameter map')

    # Options
    B1_map = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    f0_map = File(desc='f0 map (Hertz) file', argstr='--f0=%s')
    T2b = traits.Float(desc='T2 of bound pool', argstr='--T2b=%f')


class eMT(base.FitCommand):
    """
    Fit MT parameters to ellipse parameters

    """

    _cmd = 'qi ssfp_emt'
    input_spec = eMTInputSpec
    output_spec = base.FitOutputSpec(
        'EMT', ['PD', 'T1_f', 'T2_f', 'k_bf', 'f_b'])


class eMTSim(base.SimCommand):
    """
    Simulate ellipse parameters from MT parameters

    """

    _cmd = 'qi ssfp_emt'
    input_spec = base.SimInputSpec('EMT',
                                   varying=['PD', 'f_b',
                                            'k_bf', 'T1_f', 'T2_f'],
                                   fixed=['f0', 'B1'],
                                   out_files=['G', 'a', 'b'],
                                   extras={'T2b': traits.Float(desc='T2 of bound pool', argstr='--T2b=%f')})
    output_spec = base.SimOutputSpec('EMT', ['G', 'a', 'b'])

############################ qi_mtsat ############################


class MTSatInputSpec(base.InputSpec):
    # Inputs
    PDw_file = File(exists=True, argstr='%s', mandatory=True,
                    position=-3, desc='Path to PD-weighted data')

    T1w_file = File(exists=True, argstr='%s', mandatory=True,
                    position=-2, desc='Path to T1-weighted data')

    MTw_file = File(exists=True, argstr='%s', mandatory=True,
                    position=-1, desc='Path to MT-weighted data')

    # Options
    residuals = traits.Bool(
        desc='Write out residuals for each data-point', argstr='--resids')


class MTSat(base.FitCommand):
    """
    Runs qi_mtsat

    """

    _cmd = 'qi mtsat'
    input_spec = MTSatInputSpec
    output_spec = base.FitOutputSpec('MTSat', ['PD', 'R1', 'delta'])


class MTSatSim(base.SimCommand):
    _cmd = 'qi mtsat'
    input_spec = base.SimInputSpec('MTSat', varying=['PD', 'R1', 'delta'], fixed=[
                                   'B1'], out_files=['PDw', 'T1w', 'MTw'])
    output_spec = base.SimOutputSpec('MTSat', ['PDw', 'T1w', 'MTw'])

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
