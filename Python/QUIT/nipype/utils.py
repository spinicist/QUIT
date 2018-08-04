#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT utilities.

Currently implemented:
    - qisignal

To be implemented:
    - qi_coil_combine
    - qi_gradient
    - qi_rfprofile
    - qiaffine
    - qicomplex
    - qihdr
    - qikfilter
    - qimask
    - qipolyfit
    - qipolyimg
    - qisplitsubjects

Requires that the QUIT tools are in your your system path

"""

from __future__ import (print_function, division, unicode_literals,
                        absolute_import)

from nipype.interfaces.base import CommandLineInputSpec, CommandLine, TraitedSpec, File, traits, isdefined
import json

############################ qisignal ############################

class QiSignalInputSpec(CommandLineInputSpec):
    
    # Inputs
    param_file = File(desc='Signal Equation .json file', position=-1, argstr='< %s', 
        xor=['param_dict'], mandatory=True, exists=True)

    param_dict = traits.Dict(desc='Signal Equation Dict', position=-1, 
        argstr='', mandatory=True, xor=['param_file'])
        
    out_file = File(exists=False, argstr='%s', mandatory=True,
        position=-2, desc='Output file')
    
    # Options
    verbose = traits.Bool(desc='Print more information', argstr='-v')
    threads = traits.Int(desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(desc='Add a prefix to output filenames', argstr='--out=%s')
    mask_file = File(desc='Only process voxels within the mask', argstr='--mask=%s')
    ref_file = File(desc='Resample inputs to this reference space', argstr='--ref=%s')
    sim_noise = traits.Float(desc='Add complex noise with std=value', argstr='---noise=%f')
    sim_seed = traits.Int(desc='Seed noise RNG with specific value', argstr='--seed=%d')
    n_comp = traits.Int(desc='Choose number of components in model 1/2/3, default 1', argstr='--model=%d')
    complex_output = traits.Bool(desc='Save complex images', argstr='--complex')
    environ = {'QUIT_EXT':'NIFTI_GZ'}
    
class QiSignalOutputSpec(TraitedSpec):
    output_img = File(desc="Simulated Image")

class QiSignal(CommandLine):
    """
    Run signal simulation with qisignal

    Example with parameter file
    -------
    >>> from QUIT.nipype.utilities import QiSignal
    >>> qisignal = QiSignal(out_file='test.nii', param_file='seq_params.json')
    >>> sim_res = d1.run()
    >>> print(sim_res.outputs)

    Example with parameter dictionary
    -------
    >>> from QUIT.nipype.relaxometry import QiDespot1
    >>> params = {'SPGR': {'TR':5E-3, 'FA':[5,10]} }
    >>> d = {
            "PD": "pd.nii",
            "T1": "t1.nii.gz",
            "T2": "",
            "f0": "",
            "B1": "",
            "SequenceGroup": {
                "sequences": [
                    {
                        "SPGR": {
                            "TR": 0.05,
                            "FA": [3, 10]
                            }
                        }
                    ]
                }
            }
    >>> qisignal = QiSignal(out_file='test.nii', param_dict=d)
    >>> sim_res = d1.run()
    >>> print(sim_res.outputs)

    """

    _cmd = 'qisignal'
    input_spec = QiSignalInputSpec
    output_spec = QiSignalOutputSpec

    def _format_arg(self, name, spec, value):
        if name == 'param_dict':
            with open('_tmp_input.json', 'w') as outfile:
                json.dump(value, outfile)
            return "< _tmp_input.json"

        return super(QiSignal, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        
        prefix = ''
        if self.inputs.prefix:
            prefix = self.inputs.prefix
            
        outputs['output_img'] = prefix + self.inputs.out_file
        
        return outputs