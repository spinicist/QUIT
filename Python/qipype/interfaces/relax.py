#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT relaxometry utilities
Requires that the QUIT tools are in your your system path
"""

from os import path, getcwd
from nipype.interfaces.base import TraitedSpec, File, traits
from .. import base


############################ qidespot1 ############################

class DESPOT1InputSpec(base.FitInputSpec):
    # Additional Options
    b1_map = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    algo = traits.String(desc="Choose algorithm (l/w/n)", argstr="--algo=%s")
    iterations = traits.Int(
        desc='Max iterations for WLLS/NLLS (default 15)', argstr='--its=%d')


class DESPOT1(base.FitCommand):
    """
    Run DESPOT1 analysis with qidespot1

    Example with parameter file
    -------
    >>> from quit.nipype.relaxometry import QIDespot1
    >>> d1 = QIDespot1(prefix='nipype_', param_file='spgr_params.json')
    >>> d1.inputs.in_file = 'SPGR.nii.gz'
    >>> d1_res = d1.run()
    >>> print(d1_res.outputs)

    Example with parameter dictionary
    -------
    >>> from quit.nipype.relaxometry import QIDespot1
    >>> params = {'SPGR': {'TR':5E-3, 'FA':[5,10]} }
    >>> d1 = QIDespot1(prefix='nipype_', param_dict=params)
    >>> d1.inputs.in_file = 'SPGR.nii.gz'
    >>> d1_res = d1.run()
    >>> print(d1_res.outputs)

    """

    _cmd = 'qi despot1'
    input_spec = DESPOT1InputSpec
    output_spec = base.FitOutputSpec('DESPOT1OutputSpec', 'D1', ['PD', 'T1'])

############################ qidespot1sim ############################


class DESPOT1SimInputSpec(base.SimInputSpec):
    # Options
    b1_map = File(desc='B1 map (ratio) file', argstr='--B1=%s')


class DESPOT1Sim(base.SimCommand):
    """
    Run DESPOT1 simulation with qidespot1

    Example with parameter dictionary
    -------
    >>> from quit.nipype.relaxometry import QIDespot1Sim
    >>> params = {'SPGR': {'TR':5E-3, 'FA':[5,10]},
                  'T1File': 'D1_T1.nii.gz',
                  'PDFile': 'D1_PD.nii.gz'}
    >>> d1sim = QIDespot1Sim(prefix='nipype_', param_dict=params)
    >>> d1sim.inputs.spgr_file = 'SPGR.nii.gz'
    >>> d1sim_res = d1.run()
    >>> print(d1sim_res.outputs)

    """

    _cmd = 'qi despot1'
    _param_files = ['PD', 'T1']
    input_spec = DESPOT1SimInputSpec
    output_spec = base.SimOutputSpec

############################ qidespot1hifi ############################


class HIFIInputSpec(base.InputSpec):
    # Inputs
    spgr_file = File(exists=True, argstr='%s', mandatory=True,
                     position=-2, desc='Path to SPGR data')

    mprage_file = File(exists=True, argstr='%s', mandatory=True,
                       position=-1, desc='Path to MPRAGE data')

    # Options
    residuals = traits.Bool(
        desc='Write out residuals for each data-point', argstr='--resids')
    clamp_T1 = traits.Float(
        desc='Clamp T1 between 0 and value', argstr='--clampT1=%f')


class HIFIOutputSpec(TraitedSpec):
    t1_map = File('HIFI_T1.nii.gz', desc="Path to T1 map", usedefault=True)
    pd_map = File('HIFI_PD.nii.gz', desc="Path to PD map", usedefault=True)
    b1_map = File('HIFI_B1.nii.gz', desc="Path to B1 map", usedefault=True)
    rmse_map = File('HIFI_rmse.nii.gz',
                    desc="Path to residual map", usedefault=True)


class HIFI(base.FitCommand):
    """
    Calculate T1 & B1 map with the DESPOT1-HIFI method

    Example
    -------
    >>> from quit.nipype.relaxometry import QIDespot1Hifi
    >>> params = {'SPGR': {'TR':5E-3, 'FA':[5,10]}, 
                  'MPRAGE': { 'FA': 5, 'TR': 5E-3, 'TI': 0.45, 'TD': 0, 'eta': 1, 'ETL': 64, 'k0': 0 }}
    >>> hifi = QIDespot1Hifi(prefix='nipype_', param_dict=params)
    >>> hifi.inputs.spgr_file = 'SPGR.nii.gz'
    >>> hifi.inputs.mprage_file = 'MPRAGE.nii.gz'
    >>> hifi_res = hifi.run()
    >>> print(hifi_res.outputs)

    """

    _cmd = 'qi despot1hifi'
    input_spec = HIFIInputSpec
    output_spec = HIFIOutputSpec

############################ qidespot1hifisim ############################


class HIFISimInputSpec(base.SimInputBaseSpec):
    # Inputs
    spgr_file = File(exists=False, argstr='%s', mandatory=True,
                     position=-2, desc='Path for output SPGR image')

    mprage_file = File(exists=False, argstr='%s', mandatory=True,
                       position=-1, desc='Path for output MPRAGE image')


class HIFISimOutputSpec(TraitedSpec):
    spgr_file = File(desc="Output SPGR image")
    mprage_file = File(desc="Output MPRAGE image")


class HIFISim(base.SimCommand):
    """
    Simulate SPGR/FLASH and MPRAGE images using DESPOT1-HIFI model

    Example
    -------
    >>> from quit.nipype.relaxometry import QIDespot1HifiSim
    >>> params = {'SPGR': {'TR':5E-3, 'FA':[5,10]}, 
                  'MPRAGE': { 'FA': 5, 'TR': 5E-3, 'TI': 0.45, 'TD': 0, 'eta': 1, 'ETL': 64, 'k0': 0 },
                  'PDFile': 'PD.nii.gz',
                  'T1File': 'T1.nii.gz'}
    >>> hifisim = QIDespot1HifiSim(prefix='nipype_', param_dict=params)
    >>> hifisim.inputs.spgr_file = 'SPGR.nii.gz'
    >>> hifisim.inputs.mprage_file = 'MPRAGE.nii.gz'
    >>> hifisim_res = hifi.run()
    >>> print(hifisim_res.outputs)

    """

    _cmd = 'qi despot1hifi'
    _param_files = ['PD', 'T1', 'B1']
    input_spec = HIFISimInputSpec
    output_spec = HIFISimOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['spgr_file'] = self.inputs.spgr_file
        outputs['mprage_file'] = self.inputs.mprage_file
        return outputs

############################ qidespot2 ############################


class DESPOT2InputSpec(base.FitInputSpec):
    # Inputs
    t1_file = File(exists=True, argstr='%s', mandatory=True,
                   position=-2, desc='Path to T1 map')

    # Options
    b1_map = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    algo = traits.Enum("LLS", "WLS", "NLS",
                       desc="Choose algorithm", argstr="--algo=%d")
    ellipse = traits.Bool(
        desc="Data is ellipse geometric solution", argstr='--gs')
    iterations = traits.Int(
        desc='Max iterations for WLLS/NLLS (default 15)', argstr='--its=%d')
    clamp_PD = traits.Float(
        desc='Clamp PD between 0 and value', argstr='-f %f')
    clamp_T2 = traits.Float(
        desc='Clamp T2 between 0 and value', argstr='--clampT1=%f')


class DESPOT2OutputSpec(TraitedSpec):
    t2_map = File('D2_T2.nii.gz', desc="Path to T2 map", usedefault=True)
    pd_map = File('D2_PD.nii.gz', desc="Path to PD map", usedefault=True)
    rmse_map = File('D2_rmse.nii.gz',
                    desc="Path to residual map", usedefault=True)


class DESPOT2(base.FitCommand):
    """
    Run DESPOT2 analysis with qidespot2

    Example with parameter file
    -------
    >>> from quit.nipype.relaxometry import QIDespot2
    >>> d1 = QIDespot2(prefix='nipype_', param_file='ssfp_params.json')
    >>> d2.inputs.in_file = 'SSFP.nii.gz'
    >>> d2.inputs.t1_file = 'D1_T1.nii.gz'
    >>> d2_res = d2.run()
    >>> print(d2_res.outputs)

    Example with parameter dictionary
    -------
    >>> from quit.nipype.relaxometry import QIDespot2
    >>> params = {'SSFP': {'TR':5E-3, 'FA':[10,50]} }
    >>> d2 = QIDespot2(prefix='nipype_', param_dict=params)
    >>> d2.inputs.in_file = 'SSFP.nii.gz'
    >>> d2.inputs.t1_file = 'D1_T1.nii.gz'
    >>> d2_res = d2.run()
    >>> print(d2_res.outputs)

    """

    _cmd = 'qi despot2'
    input_spec = DESPOT2InputSpec
    output_spec = DESPOT2OutputSpec

############################ qidespot2sim ############################


class DESPOT2SimInputSpec(base.SimInputSpec):
    # Inputs
    t1_file = File(exists=True, argstr='%s', mandatory=True,
                   position=-2, desc='Path for input T1 map')

    # Options
    b1_map = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    ellipse = traits.Bool(
        desc="Data is ellipse geometric solution", argstr='--gs')


class DESPOT2Sim(base.SimCommand):
    """
    Run DESPOT2 simulation

    Example with parameter dictionary
    -------
    >>> from quit.nipype.relaxometry import QIDespot2Sim
    >>> params = {'SSFP': {'TR':5E-3, 'FA':[10,50]},
                  'PDFile': 'PD.nii.gz',
                  'T1File': 'T1.nii.gz' }
    >>> d2sim = QIDespot2Sim(prefix='nipype_', param_dict=params)
    >>> d2sim.inputs.in_file = 'SSFP.nii.gz'
    >>> d2sim.inputs.t1_file = 'D1_T1.nii.gz'
    >>> d2sim_res = d2.run()
    >>> print(d2sim_res.outputs)

    """

    _cmd = 'qi despot2'
    _param_files = ['PD', 'T2']
    input_spec = DESPOT2SimInputSpec
    output_spec = base.SimOutputSpec

############################ qidespot2fm ############################


class FMInputSpec(base.FitInputSpec):
    # Inputs
    t1_file = File(exists=True, argstr='%s', mandatory=True,
                   position=-2, desc='Path to T1 map')

    # Options
    b1_map = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    asym = traits.Bool(
        desc="Fit asymmetric (+/-) off-resonance frequency", argstr='--asym')
    algo = traits.Enum("LLS", "WLS", "NLS",
                       desc="Choose algorithm", argstr="--algo=%d")


class FMOutputSpec(TraitedSpec):
    t2_map = File('FM_T2.nii.gz', desc="Path to T2 map", usedefault=True)
    pd_map = File('FM_PD.nii.gz', desc="Path to PD map", usedefault=True)
    f0_map = File('FM_f0.nii.gz',
                  desc="Path to f0 (off-resonance) map", usedefault=True)
    rmse_map = File('FM_rmse.nii.gz',
                    desc="Path to residual map", usedefault=True)


class FM(base.FitCommand):
    """
    Run DESPOT2-FM analysis

    Example
    -------
    >>> from quit.nipype.relaxometry import QIDespot2FM
    >>> params = {'SSFP': {'TR':5E-3, 'FA':[10,10,50,50], 'PhaseInc':[180,180,0,0] }
    >>> fm = QIDespot2FM(prefix='nipype_', param_dict=params)
    >>> fm.inputs.in_file = 'SSFP.nii.gz'
    >>> fm.inputs.t1_file = 'D1_T1.nii.gz'
    >>> fm_res = fm.run()
    >>> print(fm_res.outputs)

    """

    _cmd = 'qi despot2fm'
    input_spec = FMInputSpec
    output_spec = FMOutputSpec


class FMSimInputSpec(base.SimInputSpec):
    # Inputs
    t1_file = File(exists=True, argstr='%s', mandatory=True,
                   position=-2, desc='Path for input T1 map')

    # Options
    b1_map = File(desc='B1 map (ratio)', argstr='--B1=%s')
    asym = traits.Bool(
        desc="Fit asymmetric (+/-) off-resonance frequency", argstr='--asym')


class FMSim(base.SimCommand):
    """
    Run DESPOT2-FM simulation

    """

    _cmd = 'qi despot2fm'
    _param_files = ['PD', 'T2', 'f0']
    input_spec = FMSimInputSpec
    output_spec = base.SimOutputSpec

############################ qi_jsr ############################


class JSRInputSpec(base.InputSpec):
    # Inputs
    spgr_file = File(exists=True, argstr='%s', mandatory=True,
                     position=-2, desc='Path to SPGR data')

    ssfp_file = File(exists=True, argstr='%s', mandatory=True,
                     position=-1, desc='Path to SSFP data')

    # Options
    b1_map = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    residuals = traits.Bool(
        desc='Write out residuals for each data-point', argstr='--resids')
    npsi = traits.Int(
        desc='Number of psi/off-resonance starts', argstr='--npsi=%d')


class JSROutputSpec(TraitedSpec):
    pd_map = File('JSR_PD.nii.gz', desc='Path to PD map', usedefault=True)
    t1_map = File('JSR_T1.nii.gz', desc='Path to T1 map', usedefault=True)
    t2_map = File('JSR_T2.nii.gz', desc='Path to T2 map', usedefault=True)
    f0_map = File('JSR_df0.nii.gz',
                  desc='Path to off-resonance map', usedefault=True)
    rmse_map = File('JSR_rmse.nii.gz',
                    desc='Path to residual map', usedefault=True)


class JSR(base.FitCommand):
    """
    Calculate T1 &T2 map with Joint System Relaxometry

    """

    _cmd = 'qi jsr'
    input_spec = JSRInputSpec
    output_spec = JSROutputSpec

############################ qimcdespot ############################
# Status: Everything is there but not tested


class QIMCDespotInputSpec(base.InputSpec):
    # Inputs
    spgr_file = File(exists=True, argstr='%s', mandatory=True,
                     position=0, desc='SPGR file')
    ssfp_file = File(exists=True, argstr='%s', mandatory=True,
                     position=1, desc='SSFP file')
    param_dict = traits.Dict(desc='Parameters as dictionary', position=2,
                             argstr='', mandatory=True, xor=['param_file'])
    param_file = File(desc='Parameters as .json file', position=2, argstr='--json=%s',
                      xor=['param_dict'], mandatory=True, exists=True)

    # Options
    f0_map = File(desc='f0 map (Hertz)', argstr='--f0=%s', exists=True)
    B1_map = File(desc='B1 map (ratio)', argstr='--B1=%s', exists=True)
    mask = File(desc='Only process voxels within the mask',
                argstr='--mask=%s', exists=True)
    residuals = traits.Bool(
        desc='Write out residuals for each data-point', argstr='--resid')
    model = traits.Enum("1", "2", "2nex", "3", "3_f0", "3nex", desc="Select model to fit - 1/2/2nex/3/3_f0/3nex, default 3",
                        argstr="--model=%d")
    scale = traits.Bool(
        desc='Normalize signals to mean (a good idea)', argstr='--scale')
    algo = traits.Enum(
        "S", "G", desc="Select (S)tochastic or (G)aussian Region Contraction", argstr="--algo=%d")
    iterations = traits.Int(
        desc='Max iterations, default 4', argstr='---its=%d')
    field_strength = traits.Float(
        desc='Specify field-strength for fitting regions - 3/7/u for user input', argstr='--tesla=%f')
    threads = traits.Int(
        desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(
        desc='Add a prefix to output filenames', argstr='--out=%s')


class QIMCDespotOutputSpec(TraitedSpec):
    # Specify which outputs there are
    output_3C_T1_m = File(desc="T1 of myelin water")
    output_3C_T2_m = File(desc="T2 of myelin water")
    output_3C_T1_ie = File(desc="T1 of intra/extra-cellular water")
    output_3C_T2_ie = File(desc="T2 of intra/extra-cellular water")
    output_3C_T1_csf = File(desc="T1 of CSF")
    output_3C_T2_csf = File(desc="T2 of CSF")
    output_3C_tau_m = File(
        desc="The residence time of myelin water (reciprocal of forward exchange rate)")
    output_3C_f_m = File(desc="The Myelin Water Fraction (MWF)")
    output_3C_f_csf = File(desc="The CSF Fraction")
    output_3C_f0 = File(
        desc="The off-resonance frequency. If this was specified on the command line, it will be a copy of that file")
    output_3C_B1 = File(
        desc="The relative flip-angle map. If this was specified on the command line, it will be a copy of that file")


class QIMCDespot(base.FitCommand):
    """
    Interace for qimcdespot

    Example 1
    -------
    >>> from quit.nipype.relaxometry import QiMcDespot
    >>> interface = QiMcDespot(prefix='nipype_', param_file='mcdespot_params.json')
    """

    _cmd = 'qi mcdespot'
    input_spec = QIMCDespotInputSpec
    output_spec = QIMCDespotOutputSpec

    # If the command requires a json input file
    def _format_arg(self, name, spec, value):
        return self._process_params(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()

        # Specify output files
        outputs = ['3C_T1_m', '3C_T2_m', '3C_T1_ie', '3C_T2_ie', '3C_T1_csf', '3C_T2_csf',
                   '3C_tau_m', '3C_f_m', '3C_f_csf', '3C_f0', '3C_B1']
        for op in outputs:
            outputs['output_{}'.format(op)] = os.path.abspath(
                self._add_prefix('{}.nii.gz'.format(op)))
        return outputs

############################ qimp2rage ############################
# Implemented but not tested #


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
    t1_map = File('MP2_T1.nii.gz', desc='T1 Map', usedefault=True)


class MP2RAGE(base.FitCommand):
    """
    Interface for qimp2rage

    Example 1
    -------
    >>> from quit.nipype.relaxometry import QIMP2RAGE
    >>> interface = QIMP2RAGE(prefix='nipype_', param_file='spgr_params.json')

    """

    _cmd = 'qi mp2rage'
    input_spec = MP2RAGEInputSpec
    output_spec = MP2RAGEOutputSpec

############################ qimultiecho ############################


class MultiechoInputSpec(base.FitInputSpec):
    # Options
    algo = traits.String(desc="Choose algorithm (l/a/n)", argstr="--algo=%s")
    iterations = traits.Int(
        desc='Max iterations for WLLS/NLLS (default 15)', argstr='--its=%d')
    thresh_PD = traits.Float(
        desc='Only output maps when PD exceeds threshold value', argstr='-t=%f')
    clamp_T2 = traits.Float(
        desc='Clamp T2 between 0 and value', argstr='-p=%f')


class MultiechoOutputSpec(TraitedSpec):
    t2_map = File(
        'ME_T2.nii.gz', desc='The T2 map. Units are the same as TE1 and ESP', usedefault=True)
    pd_map = File(
        'ME_PD.nii.gz', desc='The apparent proton-density map (intercept of the decay curve at TE=0)', usedefault=True)
    rmse_map = File('ME_rmse.nii.gz',
                    desc='Path to residual map', usedefault=True)


class Multiecho(base.FitCommand):
    """
    help for myInterface

    Example 1
    -------
    >>> from quit.nipype.relaxometry import QiMultiecho
    >>> interface = QiMultiecho(prefix='nipype_', param_file='me_params.json')

    """

    _cmd = 'qi multiecho'
    input_spec = MultiechoInputSpec
    output_spec = MultiechoOutputSpec


class MultiechoSim(base.SimCommand):
    """
    help for myInterface

    Example 1
    -------
    >>> from quit.nipype.relaxometry import QiMultiecho
    >>> interface = QiMultiecho(prefix='nipype_', param_file='me_params.json')

    """

    _cmd = 'qi multiecho'
    _param_files = ['PD', 'T2']
    input_spec = base.SimInputSpec
    output_spec = base.SimOutputSpec

############################ qi_mpm_r2s ############################


class MPMR2sInputSpec(base.InputSpec):
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


class MPMR2sOutputSpec(TraitedSpec):
    r2s_map = File(
        'MPM_R2s.nii.gz', desc='The R2* map. Units are the same as TE1 and ESP', usedefault=True)
    s0_pdw = File('MPM_S0_PDw.nii.gz',
                  desc='Intercept of the decay curve at TE=0 for PDw', usedefault=True)
    s0_t1w = File('MPM_S0_T1w.nii.gz',
                  desc='Intercept of the decay curve at TE=0 for T1w', usedefault=True)
    s0_mtw = File('MPM_S0_MTw.nii.gz',
                  desc='Intercept of the decay curve at TE=0 for MTw', usedefault=True)


class MPMR2s(base.FitCommand):
    """
    Runs qi_mpm_r2s

    """

    _cmd = 'qi mpm_r2s'
    input_spec = MPMR2sInputSpec
    output_spec = MPMR2sOutputSpec

############################ qidream ############################
# Implemented but not tested #


class QIDreamInputSpec(base.InputSpec):
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


class QIDreamOutputSpec(TraitedSpec):
    b1_rel_map = File(desc="The relative flip-angle map.")
    b1_act_map = File(desc="The actual achieved angle in each voxel.")


class QIDream(base.BaseCommand):
    """
    Interface for qidream

    Example 1
    -------
    >>> from quit.nipype.relaxometry import QiDream
    >>> interface = QiDream(prefix='nipype_', param_file='spgr_params.json')

    """

    _cmd = 'qi dream'
    input_spec = QIDreamInputSpec
    output_spec = QIDreamOutputSpec

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


class QIAFIInputSpec(base.InputSpec):
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


class QIAFIOutputSpec(TraitedSpec):
    # Specify which outputs there are
    b1_rel_map = File(desc="The relative flip-angle map.")
    b1_act_map = File(desc="The actual flip-angle map.")


class QIAFI(base.BaseCommand):
    """
    help for myInterface

    Example 1
    -------
    >>> from quit.nipype.relaxometry import QiAfi
    >>> interface = QiAfi(prefix='nipype_', param_file='spgr_params.json')

    """

    _cmd = 'qi afi'
    input_spec = QIAFIInputSpec
    output_spec = QIAFIOutputSpec

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

############################ qi_ssfp_elipse ############################


class EllipseInputSpec(base.FitInputSpec):
    # Additional Options
    algo = traits.String(desc='Choose algorithm (h/d)', argstr='--algo=%s')


class EllipseOutputSpec(TraitedSpec):
    G_map = File('ES_G.nii.gz', desc='Path to G map', usedefault=True)
    a_map = File('ES_a.nii.gz', desc='Path to a map', usedefault=True)
    b_map = File('ES_b.nii.gz', desc='Path to b map', usedefault=True)
    theta0_map = File('ES_theta_0.nii.gz',
                      desc='Path to theta 0 (off-resonance phase) map', usedefault=True)
    phi_rf_map = File('ES_phi_rf.nii.gz',
                      desc='Path to RF phase map', usedefault=True)
    rmse_map = File('ES_rmse.nii.gz',
                    desc='Path to residual map', usedefault=True)


class Ellipse(base.FitCommand):
    """
    Fit an ellipse to SSFP data

    """

    _cmd = 'qi ssfp_ellipse'
    input_spec = EllipseInputSpec
    output_spec = EllipseOutputSpec


class EllipseSim(base.SimCommand):
    """
    Simulate SSFP data from ellipse parameters

    """

    _cmd = 'qi ssfp_ellipse'
    _param_files = ['G', 'a', 'b', 'theta_0', 'phi_rf']
    input_spec = base.SimInputSpec
    output_spec = base.SimOutputSpec

############################ qi_planet ############################


class PLANETInputSpec(base.InputSpec):
    # Inputs
    G_map = File(exists=True, argstr='%s', mandatory=True,
                 position=-3, desc='Path to G parameter map')
    a_map = File(exists=True, argstr='%s', mandatory=True,
                 position=-2, desc='Path to a parameter map')
    b_map = File(exists=True, argstr='%s', mandatory=True,
                 position=-1, desc='Path to b parameter map')

    # Options
    b1_map = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    residuals = traits.Bool(
        desc='Write out residuals for each data-point', argstr='--resids')


class PLANETOutputSpec(TraitedSpec):
    PD_map = File('PLANET_PD.nii.gz', desc="Path to PD map", usedefault=True)
    T1_map = File('PLANET_T1.nii.gz', desc="Path to T1 map", usedefault=True)
    T2_map = File('PLANET_T2.nii.gz', desc="Path to T2 map", usedefault=True)
    rmse_map = File('PLANET_rmse.nii.gz',
                    desc="Path to residual map", usedefault=True)


class PLANET(base.FitCommand):
    """
    Calculate T1/T2 from ellipse parameters

    """

    _cmd = 'qi planet'
    input_spec = PLANETInputSpec
    output_spec = PLANETOutputSpec


class PLANETSimInputSpec(base.SimInputBaseSpec):
    G_file = File(argstr='%s', mandatory=True,
                  position=-3, desc='Output G file')
    a_file = File(argstr='%s', mandatory=True,
                  position=-2, desc='Output a file')
    b_file = File(argstr='%s', mandatory=True,
                  position=-1, desc='Output b file')
    b1_map = File(desc='B1 map (ratio) file', argstr='--B1=%s')


class PLANETSimOutputSpec(TraitedSpec):
    G_file = File(desc='Output G file')
    a_file = File(desc='Output a file')
    b_file = File(desc='Output b file')


class PLANETSim(base.SimCommand):
    """
    Simulate ellipse parameters from T1/T2

    """

    _cmd = 'qi planet'
    _param_files = ['PD', 'T1', 'T2']
    input_spec = PLANETSimInputSpec
    output_spec = PLANETSimOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['G_file'] = path.abspath(self.inputs.G_file)
        outputs['a_file'] = path.abspath(self.inputs.a_file)
        outputs['b_file'] = path.abspath(self.inputs.b_file)
        return outputs
