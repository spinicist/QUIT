#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT relaxometry utilities
Requires that the QUIT tools are in your your system path
"""

from nipype.interfaces.base import TraitedSpec, File, traits
from . import base as QI


############################ qidespot1 ############################

class DESPOT1InputSpec(QI.FitInputSpec):
    # Additional Options
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    algo = traits.String(desc="Choose algorithm (l/w/n)", argstr="--algo=%s")
    iterations = traits.Int(desc='Max iterations for WLLS/NLLS (default 15)', argstr='--its=%d')
    clamp_PD = traits.Float(desc='Clamp PD between 0 and value', argstr='-f %f')
    clamp_T1 = traits.Float(desc='Clamp T1 between 0 and value', argstr='--clampT1=%f')


class DESPOT1OutputSpec(TraitedSpec):
    t1_map = File('D1_T1.nii.gz', desc="Path to T1 map", usedefault=True)
    pd_map = File('D1_PD.nii.gz', desc="Path to PD map", usedefault=True)
    residual_map = File('D1_SoS_residual.nii.gz', desc="Path to residual map", usedefault=True)


class DESPOT1(QI.FitCommand):
    """
    Run DESPOT1 analysis with qidespot1

    Example with parameter file
    -------
    >>> from QUIT.nipype.relaxometry import QIDespot1
    >>> d1 = QIDespot1(prefix='nipype_', param_file='spgr_params.json')
    >>> d1.inputs.in_file = 'SPGR.nii.gz'
    >>> d1_res = d1.run()
    >>> print(d1_res.outputs)

    Example with parameter dictionary
    -------
    >>> from QUIT.nipype.relaxometry import QIDespot1
    >>> params = {'SPGR': {'TR':5E-3, 'FA':[5,10]} }
    >>> d1 = QIDespot1(prefix='nipype_', param_dict=params)
    >>> d1.inputs.in_file = 'SPGR.nii.gz'
    >>> d1_res = d1.run()
    >>> print(d1_res.outputs)

    """

    _cmd = 'qidespot1'
    input_spec = DESPOT1InputSpec
    output_spec = DESPOT1OutputSpec

############################ qidespot1sim ############################


class DESPOT1SimInputSpec(QI.SimInputSpec):
    # Options
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')


class DESPOT1Sim(QI.SimCommand):
    """
    Run DESPOT1 simulation with qidespot1

    Example with parameter dictionary
    -------
    >>> from QUIT.nipype.relaxometry import QIDespot1Sim
    >>> params = {'SPGR': {'TR':5E-3, 'FA':[5,10]},
                  'T1File': 'D1_T1.nii.gz',
                  'PDFile': 'D1_PD.nii.gz'}
    >>> d1sim = QIDespot1Sim(prefix='nipype_', param_dict=params)
    >>> d1sim.inputs.spgr_file = 'SPGR.nii.gz'
    >>> d1sim_res = d1.run()
    >>> print(d1sim_res.outputs)

    """

    _cmd = 'qidespot1'
    _param_files = ['PD', 'T1']
    input_spec = DESPOT1SimInputSpec
    output_spec = QI.SimOutputSpec

############################ qidespot1hifi ############################


class HIFIInputSpec(QI.InputSpec):
    # Inputs
    spgr_file = File(exists=True, argstr='%s', mandatory=True,
                     position=-2, desc='Path to SPGR data')

    mprage_file = File(exists=True, argstr='%s', mandatory=True,
                       position=-1, desc='Path to MPRAGE data')

    # Options
    residuals = traits.Bool(desc='Write out residuals for each data-point', argstr='--resids')
    clamp_T1 = traits.Float(desc='Clamp T1 between 0 and value', argstr='--clampT1=%f')


class HIFIOutputSpec(TraitedSpec):
    t1_map = File('HIFI_T1.nii.gz', desc="Path to T1 map", usedefault=True)
    pd_map = File('HIFI_PD.nii.gz', desc="Path to PD map", usedefault=True)
    b1_map = File('HIFI_B1.nii.gz', desc="Path to B1 map", usedefault=True)
    residual_map = File('HIFI_SoS_residual.nii.gz', desc="Path to residual map", usedefault=True)


class HIFI(QI.FitCommand):
    """
    Calculate T1 & B1 map with the DESPOT1-HIFI method

    Example
    -------
    >>> from QUIT.nipype.relaxometry import QIDespot1Hifi
    >>> params = {'SPGR': {'TR':5E-3, 'FA':[5,10]}, 
                  'MPRAGE': { 'FA': 5, 'TR': 5E-3, 'TI': 0.45, 'TD': 0, 'eta': 1, 'ETL': 64, 'k0': 0 }}
    >>> hifi = QIDespot1Hifi(prefix='nipype_', param_dict=params)
    >>> hifi.inputs.spgr_file = 'SPGR.nii.gz'
    >>> hifi.inputs.mprage_file = 'MPRAGE.nii.gz'
    >>> hifi_res = hifi.run()
    >>> print(hifi_res.outputs)

    """

    _cmd = 'qidespot1hifi'
    input_spec = HIFIInputSpec
    output_spec = HIFIOutputSpec

############################ qidespot1hifisim ############################


class HIFISimInputSpec(QI.SimInputBaseSpec):
    # Inputs
    spgr_file = File(exists=False, argstr='%s', mandatory=True,
                     position=-2, desc='Path for output SPGR image')

    mprage_file = File(exists=False, argstr='%s', mandatory=True,
                       position=-1, desc='Path for output MPRAGE image')


class HIFISimOutputSpec(TraitedSpec):
    spgr_file = File(desc="Output SPGR image")
    mprage_file = File(desc="Output MPRAGE image")


class HIFISim(QI.SimCommand):
    """
    Simulate SPGR/FLASH and MPRAGE images using DESPOT1-HIFI model

    Example
    -------
    >>> from QUIT.nipype.relaxometry import QIDespot1HifiSim
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

    _cmd = 'qidespot1hifi'
    _param_files = ['PD', 'T1', 'B1']
    input_spec = HIFISimInputSpec
    output_spec = HIFISimOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['spgr_file'] = self.inputs.spgr_file
        outputs['mprage_file'] = self.inputs.mprage_file
        return outputs

############################ qidespot2 ############################


class DESPOT2InputSpec(QI.FitInputSpec):
    # Inputs
    t1_file = File(exists=True, argstr='%s', mandatory=True,
                   position=-2, desc='Path to T1 map')

    # Options
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    algo = traits.Enum("LLS", "WLS", "NLS", desc="Choose algorithm", argstr="--algo=%d")
    ellipse = traits.Bool(desc="Data is ellipse geometric solution", argstr='--gs')
    iterations = traits.Int(desc='Max iterations for WLLS/NLLS (default 15)', argstr='--its=%d')
    clamp_PD = traits.Float(desc='Clamp PD between 0 and value', argstr='-f %f')
    clamp_T2 = traits.Float(desc='Clamp T2 between 0 and value', argstr='--clampT1=%f')


class DESPOT2OutputSpec(TraitedSpec):
    t2_map = File('D2_T2.nii.gz', desc="Path to T2 map", usedefault=True)
    pd_map = File('D2_PD.nii.gz', desc="Path to PD map", usedefault=True)
    residual_map = File('D2_SoS_residual.nii.gz', desc="Path to residual map", usedefault=True)


class DESPOT2(QI.FitCommand):
    """
    Run DESPOT2 analysis with qidespot2

    Example with parameter file
    -------
    >>> from QUIT.nipype.relaxometry import QIDespot2
    >>> d1 = QIDespot2(prefix='nipype_', param_file='ssfp_params.json')
    >>> d2.inputs.in_file = 'SSFP.nii.gz'
    >>> d2.inputs.t1_file = 'D1_T1.nii.gz'
    >>> d2_res = d2.run()
    >>> print(d2_res.outputs)

    Example with parameter dictionary
    -------
    >>> from QUIT.nipype.relaxometry import QIDespot2
    >>> params = {'SSFP': {'TR':5E-3, 'FA':[10,50]} }
    >>> d2 = QIDespot2(prefix='nipype_', param_dict=params)
    >>> d2.inputs.in_file = 'SSFP.nii.gz'
    >>> d2.inputs.t1_file = 'D1_T1.nii.gz'
    >>> d2_res = d2.run()
    >>> print(d2_res.outputs)

    """

    _cmd = 'qidespot2'
    input_spec = DESPOT2InputSpec
    output_spec = DESPOT2OutputSpec

############################ qidespot2sim ############################


class DESPOT2SimInputSpec(QI.SimInputSpec):
    # Inputs
    t1_file = File(exists=True, argstr='%s', mandatory=True,
                   position=-2, desc='Path for input T1 map')

    # Options
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    ellipse = traits.Bool(desc="Data is ellipse geometric solution", argstr='--gs')


class DESPOT2Sim(QI.SimCommand):
    """
    Run DESPOT2 simulation

    Example with parameter dictionary
    -------
    >>> from QUIT.nipype.relaxometry import QIDespot2Sim
    >>> params = {'SSFP': {'TR':5E-3, 'FA':[10,50]},
                  'PDFile': 'PD.nii.gz',
                  'T1File': 'T1.nii.gz' }
    >>> d2sim = QIDespot2Sim(prefix='nipype_', param_dict=params)
    >>> d2sim.inputs.in_file = 'SSFP.nii.gz'
    >>> d2sim.inputs.t1_file = 'D1_T1.nii.gz'
    >>> d2sim_res = d2.run()
    >>> print(d2sim_res.outputs)

    """

    _cmd = 'qidespot2'
    _param_files = ['PD', 'T2']
    input_spec = DESPOT2SimInputSpec
    output_spec = QI.SimOutputSpec

############################ qidespot2fm ############################


class FMInputSpec(QI.FitInputSpec):
    # Inputs
    t1_file = File(exists=True, argstr='%s', mandatory=True,
                   position=-2, desc='Path to T1 map')

    # Options
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    asym = traits.Bool(desc="Fit asymmetric (+/-) off-resonance frequency", argstr='--asym')
    algo = traits.Enum("LLS", "WLS", "NLS", desc="Choose algorithm", argstr="--algo=%d")


class FMOutputSpec(TraitedSpec):
    t2_map = File('FM_T2.nii.gz', desc="Path to T2 map", usedefault=True)
    pd_map = File('FM_PD.nii.gz', desc="Path to PD map", usedefault=True)
    f0_map = File('FM_f0.nii.gz', desc="Path to f0 (off-resonance) map", usedefault=True)
    residual_map = File('FM_SoS_residual.nii.gz', desc="Path to residual map", usedefault=True)


class FM(QI.FitCommand):
    """
    Run DESPOT2-FM analysis

    Example
    -------
    >>> from QUIT.nipype.relaxometry import QIDespot2FM
    >>> params = {'SSFP': {'TR':5E-3, 'FA':[10,10,50,50], 'PhaseInc':[180,180,0,0] }
    >>> fm = QIDespot2FM(prefix='nipype_', param_dict=params)
    >>> fm.inputs.in_file = 'SSFP.nii.gz'
    >>> fm.inputs.t1_file = 'D1_T1.nii.gz'
    >>> fm_res = fm.run()
    >>> print(fm_res.outputs)

    """

    _cmd = 'qidespot2fm'
    input_spec = FMInputSpec
    output_spec = FMOutputSpec


class FMSimInputSpec(QI.SimInputSpec):
    # Inputs
    t1_file = File(exists=True, argstr='%s', mandatory=True,
                   position=-2, desc='Path for input T1 map')

    # Options
    b1map_file = File(desc='B1 map (ratio)', argstr='--B1=%s')
    asym = traits.Bool(desc="Fit asymmetric (+/-) off-resonance frequency", argstr='--asym')


class FMSim(QI.SimCommand):
    """
    Run DESPOT2-FM simulation

    """

    _cmd = 'qidespot2fm'
    _param_files = ['PD', 'T2', 'f0']
    input_spec = FMSimInputSpec
    output_spec = QI.SimOutputSpec

############################ qimcdespot ############################
# Status: Everything is there but not tested


class QIMCDespotInputSpec(QI.InputSpec):
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
    mask = File(desc='Only process voxels within the mask', argstr='--mask=%s', exists=True)
    residuals = traits.Bool(desc='Write out residuals for each data-point', argstr='--resid')
    model = traits.Enum("1", "2", "2nex", "3", "3_f0", "3nex", desc="Select model to fit - 1/2/2nex/3/3_f0/3nex, default 3",
                        argstr="--model=%d")
    scale = traits.Bool(desc='Normalize signals to mean (a good idea)', argstr='--scale')
    algo = traits.Enum("S", "G", desc="Select (S)tochastic or (G)aussian Region Contraction", argstr="--algo=%d")
    iterations = traits.Int(desc='Max iterations, default 4', argstr='---its=%d')
    field_strength = traits.Float(desc='Specify field-strength for fitting regions - 3/7/u for user input', argstr='--tesla=%f')
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')


class QIMCDespotOutputSpec(TraitedSpec):
    # Specify which outputs there are
    output_3C_T1_m = File(desc="T1 of myelin water")
    output_3C_T2_m = File(desc="T2 of myelin water")
    output_3C_T1_ie = File(desc="T1 of intra/extra-cellular water")
    output_3C_T2_ie = File(desc="T2 of intra/extra-cellular water")
    output_3C_T1_csf = File(desc="T1 of CSF")
    output_3C_T2_csf = File(desc="T2 of CSF")
    output_3C_tau_m = File(desc="The residence time of myelin water (reciprocal of forward exchange rate)")
    output_3C_f_m = File(desc="The Myelin Water Fraction (MWF)")
    output_3C_f_csf = File(desc="The CSF Fraction")
    output_3C_f0 = File(desc="The off-resonance frequency. If this was specified on the command line, it will be a copy of that file")
    output_3C_B1 = File(desc="The relative flip-angle map. If this was specified on the command line, it will be a copy of that file")


class QIMCDespot(QI.FitCommand):
    """
    Interace for qimcdespot

    Example 1
    -------
    >>> from QUIT.nipype.relaxometry import QiMcDespot
    >>> interface = QiMcDespot(prefix='nipype_', param_file='mcdespot_params.json')
    """

    _cmd = 'qimcdespot'
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
            outputs['output_{}'.format(op)] = os.path.abspath(self._add_prefix('{}.nii.gz'.format(op)))
        return outputs

############################ qimp2rage ############################
# Implemented but not tested #


class MP2RAGEInputSpec(QI.InputSpec):
    # Inputs
    mprage_data = File(exists=True, argstr='%s', mandatory=True,
                       position=0, desc='Path to complex MP-RAGE data')

    # Options
    mask = File(desc='Only process voxels within the mask', argstr='--mask=%s', exists=True)
    automask = traits.Bool(desc='Create a mask from the sum of squares image', argstr='--automask')

    # Commonly used options
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')


class MP2RAGEOutputSpec(TraitedSpec):
    # Specify which outputs there are
    contrast = File(desc='The MP2 contrast image. The range of this image is -0.5 to 0.5 unless the --automask option is specified, in which case it will be shifted to 0 to 1.')
    t1map = File(desc='The T1 map. Units are the same as TR and SegTR.')


class MP2RAGE(QI.FitCommand):
    """
    Interface for qimp2rage

    Example 1
    -------
    >>> from QUIT.nipype.relaxometry import QIMP2RAGE
    >>> interface = QIMP2RAGE(prefix='nipype_', param_file='spgr_params.json')

    """

    _cmd = 'qimp2rage'
    input_spec = MP2RAGEInputSpec
    output_spec = MP2RAGEOutputSpec

    # If the command requires a json input file
    def _format_arg(self, name, spec, value):
        return self._process_params(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()

        # Specify output files
        input_bname = os.path.basename(self.inputs.mprage_data)
        input_fname = os.path.splitext(os.path.split(input_bname)[-1])

        outputs['contrast'] = os.path.abspath(self._add_prefix(input_fname + '_contrast.nii.gz'))
        outputs['t1map'] = os.path.abspath(self._add_prefix(input_fname + '_T1.nii.gz'))

        return outputs

############################ qimultiecho ############################


class MultiechoInputSpec(QI.FitInputSpec):
    # Options
    algo = traits.String(desc="Choose algorithm (l/a/n)", argstr="--algo=%s")
    iterations = traits.Int(desc='Max iterations for WLLS/NLLS (default 15)', argstr='--its=%d')
    thresh_PD = traits.Float(desc='Only output maps when PD exceeds threshold value', argstr='-t=%f')
    clamp_T2 = traits.Float(desc='Clamp T2 between 0 and value', argstr='-p=%f')


class MultiechoOutputSpec(TraitedSpec):
    t2_map = File('ME_T2.nii.gz', desc='The T2 map. Units are the same as TE1 and ESP', usedefault=True)
    pd_map = File('ME_PD.nii.gz', desc='The apparent proton-density map (intercept of the decay curve at TE=0)', usedefault=True)
    residual_map = File('ME_SoS_residual.nii.gz', desc='Path to residual map', usedefault=True)


class Multiecho(QI.FitCommand):
    """
    help for myInterface

    Example 1
    -------
    >>> from QUIT.nipype.relaxometry import QiMultiecho
    >>> interface = QiMultiecho(prefix='nipype_', param_file='me_params.json')

    """

    _cmd = 'qimultiecho'
    input_spec = MultiechoInputSpec
    output_spec = MultiechoOutputSpec


class MultiechoSim(QI.SimCommand):
    """
    help for myInterface

    Example 1
    -------
    >>> from QUIT.nipype.relaxometry import QiMultiecho
    >>> interface = QiMultiecho(prefix='nipype_', param_file='me_params.json')

    """

    _cmd = 'qimultiecho'
    _param_files = ['PD', 'T2']
    input_spec = QI.SimInputSpec
    output_spec = QI.SimOutputSpec

############################ qidream ############################
# Implemented but not tested #


class QIDreamInputSpec(QI.InputSpec):
    # Inputs
    dream_file = File(exists=True, argstr='%s', mandatory=True,
                      position=0, desc='Input file. Must have 2 volumes (FID and STE)')

    # Options
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')
    order = traits.String(desc='Volume order - f/s/v - fid/ste/vst first', argstr='--order=%s')
    mask_file = File(desc='Only process voxels within the mask', argstr='--mask=%s')
    alpha = traits.Float(desc="Nominal flip-angle (default 55)", argstr="--alpha=%f")


class QIDreamOutputSpec(TraitedSpec):
    b1_rel_map = File(desc="The relative flip-angle map.")
    b1_act_map = File(desc="The actual achieved angle in each voxel.")


class QIDream(QI.BaseCommand):
    """
    Interface for qidream

    Example 1
    -------
    >>> from QUIT.nipype.relaxometry import QiDream
    >>> interface = QiDream(prefix='nipype_', param_file='spgr_params.json')

    """

    _cmd = 'qidream'
    input_spec = QIDreamInputSpec
    output_spec = QIDreamOutputSpec

    # If the command requires a json input file
    def _format_arg(self, name, spec, value):
        return self._process_params(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['b1_rel_map'] = os.path.abspath(self._add_prefix('DREAM_B1.nii.gz'))
        outputs['b1_act_map'] = os.path.abspath(self._add_prefix('DREAM_angle.nii.gz'))

        return outputs

############################ qiafi ############################
# Implemented but not tested #


class QIAFIInputSpec(QI.InputSpec):
    # Inputs
    afi_file = File(exists=True, argstr='%s', mandatory=True, position=0, desc='Input file')

    # Options
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')
    alpha = traits.Float(desc="Specify nominal flip-angle, default 55", argstr="--flip=%f")
    tr_ratio = traits.Float(desc="Specify TR2:TR1 ratio, default 5", argstr="--ratio=%f")
    save_act_b1 = traits.Bool(desc="Write out the actual flip-angle as well as B1", argstr="--save")


class QIAFIOutputSpec(TraitedSpec):
    # Specify which outputs there are
    b1_rel_map = File(desc="The relative flip-angle map.")
    b1_act_map = File(desc="The actual flip-angle map.")


class QIAFI(QI.BaseCommand):
    """
    help for myInterface

    Example 1
    -------
    >>> from QUIT.nipype.relaxometry import QiAfi
    >>> interface = QiAfi(prefix='nipype_', param_file='spgr_params.json')

    """

    _cmd = 'qiafi'
    input_spec = QIAFIInputSpec
    output_spec = QIAFIOutputSpec

    # If the command requires a json input file
    def _format_arg(self, name, spec, value):
        return self._process_params(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['b1_rel_map'] = os.path.abspath(self._add_prefix('AFI_B1.nii.gz'))
        outputs['b1_act_map'] = os.path.abspath(self._add_prefix('AFI_angle.nii.gz'))
        return outputs
