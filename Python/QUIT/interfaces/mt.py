#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT MT module
Requires that the QUIT tools are in your your system path
"""

from os import path
from nipype.interfaces.base import TraitedSpec, DynamicTraitedSpec, File, traits, isdefined
from .. import base as QI

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

    _cmd = 'qi lineshape'
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
    rmse_map = File('LTZ_rmse.nii.gz',
                    desc='Path to residual map', usedefault=True)


class Lorentzian(QI.FitCommand):
    """
    Fit a Lorentzian function to a Z-spectrum
    """

    _cmd = 'qi lorentzian'
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
    _cmd = 'qi lorentzian'
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
    rmse_map = File('QMT_rmse.nii.gz',
                    desc="Path to residual map", usedefault=True)


class qMT(QI.FitCommand):
    """
    Fit the Ramani model to a Z-spectrum
    """

    _cmd = 'qi qmt'
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
    _cmd = 'qi qmt'
    _param_files = ['RM0a', 'f_over_R_af', 'T2_b', 'T1_a_over_T2_a', 'gM0_a']
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


class eMTInputSpec(QI.InputSpec):
    # Inputs
    G_map = File(exists=True, argstr='%s', mandatory=True,
                 position=-3, desc='Path to G parameter map')
    a_map = File(exists=True, argstr='%s', mandatory=True,
                 position=-2, desc='Path to a parameter map')
    b_map = File(exists=True, argstr='%s', mandatory=True,
                 position=-1, desc='Path to b parameter map')

    # Options
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    f0map_file = File(desc='f0 map (Hertz) file', argstr='--f0=%s')
    T2b = traits.Float(desc='T2 of bound pool', argstr='--T2b=%f')


class eMTOutputSpec(TraitedSpec):
    PD_map = File('EMT_PD.nii.gz', desc="Path to PD map", usedefault=True)
    f_b_map = File('EMT_f_b.nii.gz',
                   desc="Path to bound-pool fraction map", usedefault=True)
    k_bf_map = File('EMT_k_bf.nii.gz',
                    desc="Path to exchange rate map", usedefault=True)
    T1_f_map = File('EMT_T1_f.nii.gz',
                    desc="Path to free-pool T1 map", usedefault=True)
    T2_f_map = File('EMT_T2_f.nii.gz',
                    desc="Path to free-pool T2 map", usedefault=True)
    rmse_map = File('EMT_rmse.nii.gz',
                    desc="Path to residual map", usedefault=True)


class eMT(QI.FitCommand):
    """
    Fit MT parameters to ellipse parameters

    """

    _cmd = 'qi ssfp_emt'
    input_spec = eMTInputSpec
    output_spec = eMTOutputSpec


class eMTSimInputSpec(QI.SimInputBaseSpec):
    # Options
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    f0map_file = File(desc='f0 map (Hertz) file', argstr='--f0=%s')
    T2b = traits.Float(desc='T2 of bound pool', argstr='--T2b=%f')
    G_file = File(argstr='%s', mandatory=True,
                  position=-3, desc='Output G file')
    a_file = File(argstr='%s', mandatory=True,
                  position=-2, desc='Output a file')
    b_file = File(argstr='%s', mandatory=True,
                  position=-1, desc='Output b file')


class eMTSimOutputSpec(TraitedSpec):
    G_file = File(desc='Output G file')
    a_file = File(desc='Output a file')
    b_file = File(desc='Output b file')


class eMTSim(QI.SimCommand):
    """
    Simulate ellipse parameters from MT parameters

    """

    _cmd = 'qi ssfp_emt'
    _param_files = ['PD', 'f_b', 'k_bf', 'T1_f', 'T2_f']
    input_spec = eMTSimInputSpec
    output_spec = eMTSimOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['G_file'] = path.abspath(self.inputs.G_file)
        outputs['a_file'] = path.abspath(self.inputs.a_file)
        outputs['b_file'] = path.abspath(self.inputs.b_file)
        return outputs

############################ qi_mtsat ############################


class MTSatInputSpec(QI.InputSpec):
    # Inputs
    pdw_file = File(exists=True, argstr='%s', mandatory=True,
                    position=-3, desc='Path to PD-weighted data')

    t1w_file = File(exists=True, argstr='%s', mandatory=True,
                    position=-2, desc='Path to T1-weighted data')

    mtw_file = File(exists=True, argstr='%s', mandatory=True,
                    position=-1, desc='Path to MT-weighted data')

    # Options
    residuals = traits.Bool(
        desc='Write out residuals for each data-point', argstr='--resids')


class MTSatOutputSpec(TraitedSpec):
    s0_map = File('MTSat_S0.nii.gz', desc='S0/PD Map', usedefault=True)
    r1_map = File('MTSat_R1.nii.gz', desc='R1 map', usedefault=True)
    delta_map = File('MTSat_delta.nii.gz',
                     desc='MTSat delta map', usedefault=True)


class MTSat(QI.FitCommand):
    """
    Runs qi_mtsat

    """

    _cmd = 'qi mtsat'
    input_spec = MTSatInputSpec
    output_spec = MTSatOutputSpec

############################ qi mtr ############################


class MTRInputSpec(QI.InputSpec):
    # Options
    in_file = File(exists=True, argstr='%s', mandatory=True,
                   position=-1, desc='Input file')


class MTROutputSpec(DynamicTraitedSpec):
    pass


class MTR(QI.BaseCommand):
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
