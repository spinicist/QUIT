#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT relaxometry utilities.

Currently implemented:
    - qidespot1
    - qidespot1hifi
    - qidespot2
    - qidespot2fm

To be implemented:
    - qimcdespot
    - qimp2rage
    - qimultiecho
    - qidream
    - qiafi

Requires that the QUIT tools are in your your system path

"""

from __future__ import (print_function, division, unicode_literals,
                        absolute_import)
import os

from nipype.interfaces.base import CommandLineInputSpec, CommandLine, TraitedSpec, File, traits, isdefined
import json

############################ qidespot1 ############################

class QiDespot1InputSpec(CommandLineInputSpec):
    # Inputs
    spgr_file = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Path to SPGR data')

    param_file = File(desc='Parameter .json file', position=1, argstr='< %s', 
        xor=['param_dict'], mandatory=True, exists=True)

    param_dict = traits.Dict(desc='dictionary trait', position=1, 
        argstr='', mandatory=True, xor=['param_file'])

    # Options
    verbose = traits.Bool(desc='Print more information', argstr='-v')
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    mask_file = File(desc='Only process voxels within the mask', argstr='--mask=%s')
    residuals = traits.Bool(desc='Write out residuals for each data-point', argstr='--resids')
    algo = traits.Enum("LLS","WLS","NLS", desc="Choose algorithm", argstr="--algo=%d")
    iterations = traits.Int(desc='Max iterations for WLLS/NLLS (default 15)', argstr='--its=%d')
    clamp_PD = traits.Float(desc='Clamp PD between 0 and value', argstr='-f %f')
    clamp_T1 = traits.Float(desc='Clamp T1 between 0 and value', argstr='--clampT1=%f')
    environ = {'QUIT_EXT':'NIFTI_GZ'}
    
class QiDespot1OutputSpec(TraitedSpec):
    t1_map = File(desc="Path to T1 map")
    pd_map = File(desc="Path to PD map")
    residual_map = File(desc="Path to residual map")

class QiDespot1(CommandLine):
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
        if name == 'param_dict':
            with open('_tmp_input.json', 'w') as outfile:
                json.dump(value, outfile)
            return "< _tmp_input.json"

        return super(QiDespot1, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        
        prefix = ''
        if self.inputs.prefix:
            prefix = self.inputs.prefix
            
        outputs['t1_map'] = prefix + 'D1_T1.nii.gz'
        outputs['pd_map'] = prefix + 'D1_PD.nii.gz'
        
        if self.inputs.residuals:
            outputs['residual_map'] = prefix + 'D1_residual.nii.gz'
        
        return outputs

############################ qidespot1hifi ############################

class QiDespot1HifiInputSpec(CommandLineInputSpec):
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
    verbose = traits.Bool(desc='Print more information', argstr='-v')
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    mask_file = File(desc='Only process voxels within the mask', argstr='--mask=%s')
    residuals = traits.Bool(desc='Write out residuals for each data-point', argstr='--resids')
    clamp_T1 = traits.Float(desc='Clamp T1 between 0 and value', argstr='--clampT1=%f')
    environ = {'QUIT_EXT':'NIFTI_GZ'}

class QiDespot1HifiOutputSpec(TraitedSpec):
    t1_map = File(desc="Path to T1 map")
    pd_map = File(desc="Path to PD map")
    b1_map = File(desc="Path to B1 map")
    residual_map = File(desc="Path to residual map")

class QiDespot1Hifi(CommandLine):
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
        if name == 'param_dict':
            with open('_tmp_input.json', 'w') as outfile:
                json.dump(value, outfile)
            return "< _tmp_input.json"

        return super(QiDespot1Hifi, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        
        prefix = ''
        if self.inputs.prefix:
            prefix = self.inputs.prefix
            
        outputs['t1_map'] = prefix + 'HIFI_T1.nii.gz'
        outputs['pd_map'] = prefix + 'HIFI_PD.nii.gz'
        outputs['b1_map'] = prefix + 'HIFI_B1.nii.gz'

        if self.inputs.residuals:
            outputs['residual_map'] = prefix + 'HIFI_residual.nii.gz'
        
        return outputs

############################ qidespot2 ############################

class QiDespot2InputSpec(CommandLineInputSpec):
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
    verbose = traits.Bool(desc='Print more information', argstr='-v')
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
    environ = {'QUIT_EXT':'NIFTI_GZ'}

class QiDespot2OutputSpec(TraitedSpec):
    t2_map = File(desc="Path to T2 map")
    pd_map = File(desc="Path to PD map")
    residual_map = File(desc="Path to residual map")

class QiDespot2(CommandLine):
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
        if name == 'param_dict':
            with open('_tmp_input.json', 'w') as outfile:
                json.dump(value, outfile)
            return "< _tmp_input.json"

        return super(QiDespot2, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        
        prefix = ''
        if self.inputs.prefix:
            prefix = self.inputs.prefix
            
        outputs['t2_map'] = prefix + 'D2_T2.nii.gz'
        outputs['pd_map'] = prefix + 'D2_PD.nii.gz'
        
        if self.inputs.residuals:
            outputs['residual_map'] = prefix + 'D2_residual.nii.gz'
        
        return outputs

############################ qidespot2fm ############################

class QiDespot2FMInputSpec(CommandLineInputSpec):
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
    verbose = traits.Bool(desc='Print more information', argstr='-v')
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')
    b1map_file = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    mask_file = File(desc='Only process voxels within the mask', argstr='--mask=%s')
    asym = traits.Bool(desc="Fit asymmetric (+/-) off-resonance frequency", argstr='-A')
    residuals = traits.Bool(desc='Write out residuals for each data-point', argstr='--resids')
    algo = traits.Enum("LLS","WLS","NLS", desc="Choose algorithm", argstr="--algo=%d")
    environ = {'QUIT_EXT':'NIFTI_GZ'}

class QiDespot2FMOutputSpec(TraitedSpec):
    t2_map = File(desc="Path to T2 map")
    pd_map = File(desc="Path to PD map")
    f0_map = File(desc="Path to f0 (off-resonance) map")
    residual_map = File(desc="Path to residual map")

class QiDespot2FM(CommandLine):
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
        if name == 'param_dict':
            with open('_tmp_input.json', 'w') as outfile:
                json.dump(value, outfile)
            return "< _tmp_input.json"

        return super(QiDespot2FM, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        
        prefix = ''
        if self.inputs.prefix:
            prefix = self.inputs.prefix
            
        outputs['t2_map'] = prefix + 'FM_T2.nii.gz'
        outputs['pd_map'] = prefix + 'FM_PD.nii.gz'
        outputs['f0_map'] = prefix + 'FM_f0.nii.gz'
        
        if self.inputs.residuals:
            outputs['residual_map'] = prefix + 'FM_residual.nii.gz'
        
        return outputs

############################ qimcdespot ############################
# < To be implemented > #
class QiMcDespotInputSpec(CommandLineInputSpec):
    # Inputs
    in_file = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Input file')
    # <Add more inputs here>

    # Options
    # <Add in options here>
    
    # Commonly used options
    environ = {'QUIT_EXT':'NIFTI_GZ'}
    verbose = traits.Bool(desc='Print more information', argstr='-v')

class QiMcDespotOutputSpec(TraitedSpec):
    # Specify which outputs there are
    out_file = File(desc="Output file")
    
class QiMcDespot(CommandLine):
    """
    Interace for qimcdespot

    Example 1
    -------
    >>> from QUIT.nipype.relaxometry import QiMcDespot
    >>> interface = QiMcDespot(prefix='nipype_', param_file='spgr_params.json')
    
    """

    _cmd = 'qimcdespot'
    input_spec = QiMcDespotInputSpec
    output_spec = QiMcDespotOutputSpec

    # If the command requires a json input file
    def _format_arg(self, name, spec, value):
        if name == 'param_dict':
            with open('_tmp_input.json', 'w') as outfile:
                json.dump(value, outfile)
            return "< _tmp_input.json"

        return super(QiMcDespot, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        
        prefix = ''
        if self.inputs.prefix:
            prefix = self.inputs.prefix
        
        # Specify output files
        outputs['out_file'] = os.path.abspath(prefix + 'out_file.nii.gz')
        
        return outputs

############################ qimp2rage ############################
# < To be implemented > #
class QiMp2rageInputSpec(CommandLineInputSpec):
    # Inputs
    in_file = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Input file')
    # <Add more inputs here>

    # Options
    # <Add in options here>
    
    # Commonly used options
    environ = {'QUIT_EXT':'NIFTI_GZ'}
    verbose = traits.Bool(desc='Print more information', argstr='-v')

class QiMp2rageOutputSpec(TraitedSpec):
    # Specify which outputs there are
    out_file = File(desc="Output file")
    
class QiMp2rage(CommandLine):
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
        if name == 'param_dict':
            with open('_tmp_input.json', 'w') as outfile:
                json.dump(value, outfile)
            return "< _tmp_input.json"

        return super(QiMp2rage, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        
        prefix = ''
        if self.inputs.prefix:
            prefix = self.inputs.prefix
        
        # Specify output files
        outputs['out_file'] = os.path.abspath(prefix + 'out_file.nii.gz')
        
        return outputs

############################ qimultiecho ############################
# < To be implemented > #
class QiMultiechoInputSpec(CommandLineInputSpec):
    # Inputs
    in_file = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Input file')
    # <Add more inputs here>

    # Options
    # <Add in options here>
    
    # Commonly used options
    environ = {'QUIT_EXT':'NIFTI_GZ'}
    verbose = traits.Bool(desc='Print more information', argstr='-v')

class QiMultiechoOutputSpec(TraitedSpec):
    # Specify which outputs there are
    out_file = File(desc="Output file")
    
class QiMultiecho(CommandLine):
    """
    help for myInterface

    Example 1
    -------
    >>> from QUIT.nipype.relaxometry import QiMultiecho
    >>> interface = QiMultiecho(prefix='nipype_', param_file='spgr_params.json')
    
    """

    _cmd = 'qimultiecho'
    input_spec = QiMultiechoInputSpec
    output_spec = QiMultiechoOutputSpec

    # If the command requires a json input file
    def _format_arg(self, name, spec, value):
        if name == 'param_dict':
            with open('_tmp_input.json', 'w') as outfile:
                json.dump(value, outfile)
            return "< _tmp_input.json"

        return super(QiMultiecho, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        
        prefix = ''
        if self.inputs.prefix:
            prefix = self.inputs.prefix
        
        # Specify output files
        outputs['out_file'] = os.path.abspath(prefix + 'out_file.nii.gz')
        
        return outputs

############################ qidream ############################
# < To be implemented > #
class QiDreamInputSpec(CommandLineInputSpec):
    # Inputs
    in_file = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Input file')
    # <Add more inputs here>

    # Options
    # <Add in options here>
    
    # Commonly used options
    environ = {'QUIT_EXT':'NIFTI_GZ'}
    verbose = traits.Bool(desc='Print more information', argstr='-v')

class QiDreamOutputSpec(TraitedSpec):
    # Specify which outputs there are
    out_file = File(desc="Output file")
    
class QiDream(CommandLine):
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
        if name == 'param_dict':
            with open('_tmp_input.json', 'w') as outfile:
                json.dump(value, outfile)
            return "< _tmp_input.json"

        return super(QiDream, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        
        prefix = ''
        if self.inputs.prefix:
            prefix = self.inputs.prefix
        
        # Specify output files
        outputs['out_file'] = os.path.abspath(prefix + 'out_file.nii.gz')
        
        return outputs
############################ qiafi ############################
# < To be implemented > #
class QiAfiInputSpec(CommandLineInputSpec):
    # Inputs
    in_file = File(exists=True, argstr='%s', mandatory=True,
        position=0, desc='Input file')
    # <Add more inputs here>

    # Options
    # <Add in options here>
    
    # Commonly used options
    environ = {'QUIT_EXT':'NIFTI_GZ'}
    verbose = traits.Bool(desc='Print more information', argstr='-v')

class QiAfiOutputSpec(TraitedSpec):
    # Specify which outputs there are
    out_file = File(desc="Output file")
    
class QiAfi(CommandLine):
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
        if name == 'param_dict':
            with open('_tmp_input.json', 'w') as outfile:
                json.dump(value, outfile)
            return "< _tmp_input.json"

        return super(QiAfi, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        
        prefix = ''
        if self.inputs.prefix:
            prefix = self.inputs.prefix
        
        # Specify output files
        outputs['out_file'] = os.path.abspath(prefix + 'out_file.nii.gz')
        
        return outputs