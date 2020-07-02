#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT relaxometry tools
Requires that the QUIT tools are in your your system path
"""

from os import path, getcwd
from nipype.interfaces.base import TraitedSpec, File, traits
from .. import base


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


class MP2RAGE(base.FitCommand):
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

############################ qimcdespot ############################
# Status: Very much unsupported.


class QIMCDespot(base.FitCommand):
    """
    Interace for qimcdespot

    Example 1
    -------
    >>> from qipype.interfaces.relax import QiMcDespot
    >>> interface = QiMcDespot(prefix='nipype_', param_file='mcdespot_params.json')
    """

    _cmd = 'qi mcdespot'
    v = ['T1_m', 'T2_m', 'T1_ie', 'T2_ie', 'T1_csf',
         'T2_csf', 'tau_m', 'f_m', 'f_csf', 'f0', 'B1']
    input_spec = base.FitInputSpec('3C', fixed=['f0', 'B1'], in_files=['spgr', 'ssfp'],
                                   extra={'model': traits.Enum("1", "2", "2nex", "3", "3_f0", "3nex", desc="Select model to fit - 1/2/2nex/3/3_f0/3nex, default 3", argstr="--model=%d"),
                                          'scale': traits.Bool(desc='Normalize signals to mean (a good idea)', argstr='--scale'),
                                          'algo': traits.Enum("S", "G", desc="Select (S)tochastic or (G)aussian Region Contraction", argstr="--algo=%d"),
                                          'iterations': traits.Int(desc='Max iterations, default 4', argstr='---its=%d'),
                                          'field_strength': traits.Float(desc='Specify field-strength for fitting regions - 3/7/u for user input', argstr='--tesla=%f')}
                                   )
    output_spec = base.FitOutputSpec(
        '3C', v)
