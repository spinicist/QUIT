#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT relaxometry utilities.
Requires that the QUIT tools are in your your system path
"""

from __future__ import (print_function, division, unicode_literals,
                        absolute_import)
import os
import json

from nipype.interfaces.base import CommandLineInputSpec, CommandLine, TraitedSpec, File, traits, isdefined
from .base import QUITCommand, QUITCommandInputSpec


############################ qidespot1 ############################

class QiDespot1InputSpec(QUITCommandInputSpec):
    # Inputs
    spgr_file = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Path to SPGR data')

    param_file = File(desc='Parameter .json file', position=1, argstr='< %s', 
        xor=['param_dict'], mandatory=True, exists=True)

    param_dict = traits.Dict(desc='dictionary trait', position=1, 
        argstr='', mandatory=True, xor=['param_file'])

    # Options
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    mask_file = File(desc='Only process voxels within the mask', argstr='--mask=%s')
    residuals = traits.Bool(desc='Write out residuals for each data-point', argstr='--resids')
    algo = traits.Enum("LLS","WLS","NLS", desc="Choose algorithm", argstr="--algo=%d")
    iterations = traits.Int(desc='Max iterations for WLLS/NLLS (default 15)', argstr='--its=%d')
    clamp_PD = traits.Float(desc='Clamp PD between 0 and value', argstr='-f %f')
    clamp_T1 = traits.Float(desc='Clamp T1 between 0 and value', argstr='--clampT1=%f')
    
class QiDespot1OutputSpec(TraitedSpec):
    t1_map = File(desc="Path to T1 map")
    pd_map = File(desc="Path to PD map")
    residual_map = File(desc="Path to residual map")

class QiDespot1(QUITCommand):
    """
    Run DESPOT1 analysis with qidespot1

    Example with parameter file
    -------
    >>> from QUIT.nipype.relaxometry import QiDespot1
    >>> d1 = QiDespot1(prefix='nipype_', param_file='spgr_params.json')
    >>> d1.inputs.in_file = 'SPGR.nii.gz'
    >>> d1_res = d1.run()
    >>> print(d1_res.outputs)

    Example with parameter dictionary
    -------
    >>> from QUIT.nipype.relaxometry import QiDespot1
    >>> params = {'SPGR': {'TR':5E-3, 'FA':[5,10]} }
    >>> d1 = QiDespot1(prefix='nipype_', param_dict=params)
    >>> d1.inputs.in_file = 'SPGR.nii.gz'
    >>> d1_res = d1.run()
    >>> print(d1_res.outputs)

    """

    _cmd = 'qidespot1'
    input_spec = QiDespot1InputSpec
    output_spec = QiDespot1OutputSpec

    def _format_arg(self, name, spec, value):
        return self._process_params(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()

        outputs['t1_map'] = os.path.abspath(self._add_prefix('D1_T1.nii.gz'))
        outputs['pd_map'] = os.path.abspath(self._add_prefix('D1_PD.nii.gz'))
        
        if self.inputs.residuals:
            outputs['residual_map'] = os.path.abspath(self._add_prefix('D1_residual.nii.gz'))
        
        return outputs

############################ qidespot1hifi ############################

class QiDespot1HifiInputSpec(QUITCommandInputSpec):
    # Inputs
    spgr_file = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Path to SPGR data')

    mprage_file = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Path to MPRAGE data')

    param_file = File(desc='Parameter .json file', position=2, argstr='< %s', 
        xor=['param_dict'], mandatory=True, exists=True)

    param_dict = traits.Dict(desc='dictionary trait', position=2, 
        argstr='', mandatory=True, xor=['param_file'])

    # Options
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    mask_file = File(desc='Only process voxels within the mask', argstr='--mask=%s')
    residuals = traits.Bool(desc='Write out residuals for each data-point', argstr='--resids')
    clamp_T1 = traits.Float(desc='Clamp T1 between 0 and value', argstr='--clampT1=%f')

class QiDespot1HifiOutputSpec(TraitedSpec):
    t1_map = File(desc="Path to T1 map")
    pd_map = File(desc="Path to PD map")
    b1_map = File(desc="Path to B1 map")
    residual_map = File(desc="Path to residual map")

class QiDespot1Hifi(QUITCommand):
    """
    Calculate T1 & B1 map with the DESPOT1-HIFI method

    Example
    -------
    >>> from QUIT.nipype.relaxometry import QiDespot1Hifi
    >>> params = {'SPGR': {'TR':5E-3, 'FA':[5,10]}, 
                  'MPRAGE': { 'FA': 5, 'TR': 5E-3, 'TI': 0.45, 'TD': 0, 'eta': 1, 'ETL': 64, 'k0': 0 }}
    >>> hifi = QiDespot1Hifi(prefix='nipype_', param_dict=params)
    >>> hifi.inputs.spgr_file = 'SPGR.nii.gz'
    >>> hifi.inputs.mprage_file = 'MPRAGE.nii.gz'
    >>> hifi_res = hifi.run()
    >>> print(hifi_res.outputs)

    """

    _cmd = 'qidespot1hifi'
    input_spec = QiDespot1HifiInputSpec
    output_spec = QiDespot1HifiOutputSpec

    def _format_arg(self, name, spec, value):
        return self._process_params(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
            
        outputs['t1_map'] = os.path.abspath(self._add_prefix('HIFI_T1.nii.gz'))
        outputs['pd_map'] = os.path.abspath(self._add_prefix('HIFI_PD.nii.gz'))
        outputs['b1_map'] = os.path.abspath(self._add_prefix('HIFI_B1.nii.gz'))

        if self.inputs.residuals:
            outputs['residual_map'] = os.path.abspath(self._add_prefix('HIFI_residual.nii.gz'))
        
        return outputs

############################ qidespot2 ############################

class QiDespot2InputSpec(QUITCommandInputSpec):
    # Inputs
    t1_map = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Path to T1 map')
    ssfp_file = File(exists=True, argstr='%s', mandatory=True,
        position=1, desc='Path to SSFP data')
    param_file = File(desc='Parameter .json file', position=2, argstr='< %s', 
        xor=['param_dict'], mandatory=True, exists=True)
    param_dict = traits.Dict(desc='dictionary trait', position=2, 
        argstr='', mandatory=True, xor=['param_file'])

    # Options
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    mask_file = File(desc='Only process voxels within the mask', argstr='--mask=%s')
    residuals = traits.Bool(desc='Write out residuals for each data-point', argstr='--resids')
    algo = traits.Enum("LLS","WLS","NLS", desc="Choose algorithm", argstr="--algo=%d")
    ellipse = traits.Bool(desc="Data is ellipse geometric solution", argstr='-e')
    iterations = traits.Int(desc='Max iterations for WLLS/NLLS (default 15)', argstr='--its=%d')
    clamp_PD = traits.Float(desc='Clamp PD between 0 and value', argstr='-f %f')
    clamp_T2 = traits.Float(desc='Clamp T2 between 0 and value', argstr='--clampT1=%f')

class QiDespot2OutputSpec(TraitedSpec):
    t2_map = File(desc="Path to T2 map")
    pd_map = File(desc="Path to PD map")
    residual_map = File(desc="Path to residual map")

class QiDespot2(QUITCommand):
    """
    Run DESPOT2 analysis with qidespot2

    Example with parameter file
    -------
    >>> from QUIT.nipype.relaxometry import QiDespot2
    >>> d1 = QiDespot2(prefix='nipype_', param_file='ssfp_params.json')
    >>> d2.inputs.in_file = 'SSFP.nii.gz'
    >>> d2.inputs.t1_file = 'D1_T1.nii.gz'
    >>> d2_res = d2.run()
    >>> print(d2_res.outputs)

    Example with parameter dictionary
    -------
    >>> from QUIT.nipype.relaxometry import QiDespot2
    >>> params = {'SSFP': {'TR':5E-3, 'FA':[10,50]} }
    >>> d2 = QiDespot2(prefix='nipype_', param_dict=params)
    >>> d2.inputs.in_file = 'SSFP.nii.gz'
    >>> d2.inputs.t1_file = 'D1_T1.nii.gz'
    >>> d2_res = d2.run()
    >>> print(d2_res.outputs)

    """

    _cmd = 'qidespot2'
    input_spec = QiDespot2InputSpec
    output_spec = QiDespot2OutputSpec

    def _format_arg(self, name, spec, value):
        return self._process_params(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
            
        outputs['t2_map'] = os.path.abspath(self._add_prefix('D2_T2.nii.gz'))
        outputs['pd_map'] = os.path.abspath(self._add_prefix('D2_PD.nii.gz'))
        
        if self.inputs.residuals:
            outputs['residual_map'] = os.path.abspath(self._add_prefix('D2_residual.nii.gz'))
        
        return outputs

############################ qidespot2fm ############################

class QiDespot2FMInputSpec(QUITCommandInputSpec):
    # Inputs
    t1_map = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Path to T1 map')
    ssfp_file = File(exists=True, argstr='%s', mandatory=True,
        position=1, desc='Path to SSFP data')
    param_file = File(desc='Parameter .json file', position=2, argstr='< %s', 
        xor=['param_dict'], mandatory=True, exists=True)
    param_dict = traits.Dict(desc='dictionary trait', position=2, 
        argstr='', mandatory=True, xor=['param_file'])

    # Options
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    mask_file = File(desc='Only process voxels within the mask', argstr='--mask=%s')
    asym = traits.Bool(desc="Fit asymmetric (+/-) off-resonance frequency", argstr='-A')
    residuals = traits.Bool(desc='Write out residuals for each data-point', argstr='--resids')
    algo = traits.Enum("LLS","WLS","NLS", desc="Choose algorithm", argstr="--algo=%d")

class QiDespot2FMOutputSpec(TraitedSpec):
    t2_map = File(desc="Path to T2 map")
    pd_map = File(desc="Path to PD map")
    f0_map = File(desc="Path to f0 (off-resonance) map")
    residual_map = File(desc="Path to residual map")

class QiDespot2FM(QUITCommand):
    """
    Run DESPOT2-FM analysis

    Example
    -------
    >>> from QUIT.nipype.relaxometry import QiDespot2FM
    >>> params = {'SSFP': {'TR':5E-3, 'FA':[10,10,50,50], 'PhaseInc':[180,180,0,0] }
    >>> fm = QiDespot2FM(prefix='nipype_', param_dict=params)
    >>> fm.inputs.in_file = 'SSFP.nii.gz'
    >>> fm.inputs.t1_file = 'D1_T1.nii.gz'
    >>> fm_res = fm.run()
    >>> print(fm_res.outputs)

    """

    _cmd = 'qidespot2fm'
    input_spec = QiDespot2FMInputSpec
    output_spec = QiDespot2FMOutputSpec

    def _format_arg(self, name, spec, value):
        return self._process_params(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
    
        outputs['t2_map'] = os.path.abspath(self._add_prefix('FM_T2.nii.gz'))
        outputs['pd_map'] = os.path.abspath(self._add_prefix('FM_PD.nii.gz'))
        outputs['f0_map'] = os.path.abspath(self._add_prefix('FM_f0.nii.gz'))
        
        if self.inputs.residuals:
            outputs['residual_map'] = os.path.abspath(self._add_prefix('FM_residual.nii.gz'))
        
        return outputs

############################ qimcdespot ############################
# Status: Everything is there but not tested
class QiMcDespotInputSpec(QUITCommandInputSpec):
    # Inputs
    spgr_file = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='SPGR file')
    ssfp_file = File(exists=True, argstr='%s', mandatory=True,
        position=1, desc='SSFP file')
    param_dict = traits.Dict(desc='Parameters as dictionary', position=2, 
        argstr='', mandatory=True, xor=['param_file'])
    param_file = File(desc='Parameters as .json file', position=2, argstr='< %s', 
        xor=['param_dict'], mandatory=True, exists=True)

    # Options
    f0_map = File(desc='f0 map (Hertz)', argstr='--f0=%s', exists=True)      
    B1_map = File(desc='B1 map (ratio)', argstr='--B1=%s', exists=True)      
    mask = File(desc='Only process voxels within the mask', argstr='--mask=%s', exists=True)
    residuals = traits.Bool(desc='Write out residuals for each data-point', argstr='--resid')
    model = traits.Enum("1","2","2nex","3","3_f0","3nex", desc="Select model to fit - 1/2/2nex/3/3_f0/3nex, default 3", 
        argstr="--model=%d")
    scale = traits.Bool(desc='Normalize signals to mean (a good idea)', argstr='--scale')
    algo = traits.Enum("S","G", desc="Select (S)tochastic or (G)aussian Region Contraction", argstr="--algo=%d")
    iterations = traits.Int(desc='Max iterations, default 4', argstr='---its=%d')
    field_strength = traits.Float(desc='Specify field-strength for fitting regions - 3/7/u for user input', argstr='--tesla=%f')
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')

class QiMcDespotOutputSpec(TraitedSpec):
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
    
class QiMcDespot(QUITCommand):
    """
    Interace for qimcdespot

    Example 1
    -------
    >>> from QUIT.nipype.relaxometry import QiMcDespot
    >>> interface = QiMcDespot(prefix='nipype_', param_file='mcdespot_params.json')
    """

    _cmd = 'qimcdespot'
    input_spec = QiMcDespotInputSpec
    output_spec = QiMcDespotOutputSpec

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
class QiMp2rageInputSpec(QUITCommandInputSpec):
    # Inputs
    mprage_data = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Path to complex MP-RAGE data')
    
    # Options
    mask = File(desc='Only process voxels within the mask', argstr='--mask=%s', exists=True)
    automask = traits.Bool(desc='Create a mask from the sum of squares image', argstr='--automask')

    # Commonly used options
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')

class QiMp2rageOutputSpec(TraitedSpec):
    # Specify which outputs there are
    contrast = File(desc='The MP2 contrast image. The range of this image is -0.5 to 0.5 unless the --automask option is specified, in which case it will be shifted to 0 to 1.')
    t1map = File(desc='The T1 map. Units are the same as TR and SegTR.')
    
class QiMp2rage(QUITCommand):
    """
    Interface for qimp2rage

    Example 1
    -------
    >>> from QUIT.nipype.relaxometry import QiMp2rage
    >>> interface = QiMp2rage(prefix='nipype_', param_file='spgr_params.json')
    
    """

    _cmd = 'qimp2rage'
    input_spec = QiMp2rageInputSpec
    output_spec = QiMp2rageOutputSpec

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
# Implemented but not tested #
# Still need to put in a nice example

class QiMultiechoInputSpec(QUITCommandInputSpec):
    # Inputs
    me_file = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Input multi-echo data')

    # Options
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')
    mask_file = File(desc='Only process voxels within the mask', argstr='--mask=%s')
    algo = traits.String(desc="Choose algorithm (l/a/n)", argstr="--algo=%s")
    iterations = traits.Int(desc='Max iterations for WLLS/NLLS (default 15)', argstr='--its=%d')
    thresh_PD = traits.Float(desc='Only output maps when PD exceeds threshold value', argstr='-t=%f')
    clamp_T2 = traits.Float(desc='Clamp T2 between 0 and value', argstr='-p=%f')    

class QiMultiechoOutputSpec(TraitedSpec):
    t2_map = File(desc="The T2 map. Units are the same as TE1 and ESP.")
    pd_map = File(desc="The apparent proton-density map (intercept of the decay curve at TE=0)")
    
class QiMultiecho(QUITCommand):
    """
    help for myInterface

    Example 1
    -------
    >>> from QUIT.nipype.relaxometry import QiMultiecho
    >>> interface = QiMultiecho(prefix='nipype_', param_file='me_params.json')
    
    """

    _cmd = 'qimultiecho'
    input_spec = QiMultiechoInputSpec
    output_spec = QiMultiechoOutputSpec

    # If the command requires a json input file
    def _format_arg(self, name, spec, value):
        return self._process_params(name, spec, value)


    def _list_outputs(self):
        outputs = self.output_spec().get()
        
        outputs['t2_map'] = os.path.abspath(self._add_prefix("ME_T2.nii.gz"))
        outputs['pd_map'] = os.path.abspath(self._add_prefix("ME_PD.nii.gz"))

        return outputs

############################ qidream ############################
# Implemented but not tested #

class QiDreamInputSpec(QUITCommandInputSpec):
    # Inputs
    dream_file = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Input file. Must have 2 volumes (FID and STE)')
    
    # Options
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')
    order = traits.String(desc='Volume order - f/s/v - fid/ste/vst first', argstr='--order=%s')
    mask_file = File(desc='Only process voxels within the mask', argstr='--mask=%s')
    alpha = traits.Float(desc="Nominal flip-angle (default 55)", argstr="--alpha=%f")
    
class QiDreamOutputSpec(TraitedSpec):
    b1_rel_map = File(desc="The relative flip-angle map.")
    b1_act_map = File(desc="The actual achieved angle in each voxel.")
    
class QiDream(QUITCommand):
    """
    Interface for qidream

    Example 1
    -------
    >>> from QUIT.nipype.relaxometry import QiDream
    >>> interface = QiDream(prefix='nipype_', param_file='spgr_params.json')
    
    """

    _cmd = 'qidream'
    input_spec = QiDreamInputSpec
    output_spec = QiDreamOutputSpec

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

class QiAfiInputSpec(QUITCommandInputSpec):
    # Inputs
    afi_file = File(exists=True, argstr='%s', mandatory=True, position=0, desc='Input file')

    # Options
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')
    alpha = traits.Float(desc="Specify nominal flip-angle, default 55", argstr="--flip=%f")
    tr_ratio = traits.Float(desc="Specify TR2:TR1 ratio, default 5", argstr="--ratio=%f")
    save_act_b1 = traits.Bool(desc="Write out the actual flip-angle as well as B1", argstr="--save")

class QiAfiOutputSpec(TraitedSpec):
    # Specify which outputs there are
    b1_rel_map = File(desc="The relative flip-angle map.")
    b1_act_map = File(desc="The actual flip-angle map.")
   
class QiAfi(QUITCommand):
    """
    help for myInterface

    Example 1
    -------
    >>> from QUIT.nipype.relaxometry import QiAfi
    >>> interface = QiAfi(prefix='nipype_', param_file='spgr_params.json')
    
    """

    _cmd = 'qiafi'
    input_spec = QiAfiInputSpec
    output_spec = QiAfiOutputSpec

    # If the command requires a json input file
    def _format_arg(self, name, spec, value):
        return self._process_params(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['b1_rel_map'] = os.path.abspath(self._add_prefix('AFI_B1.nii.gz'))
        outputs['b1_act_map'] = os.path.abspath(self._add_prefix('AFI_angle.nii.gz'))
        return outputs