#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT relaxometry tools
Requires that the QUIT tools are in your your system path
"""

from os import path, getcwd
from nipype.interfaces.base import TraitedSpec, File, traits
from .. import base

DESPOT1, DESPOT1Sim = base.Command('DESPOT1', 'qi despot1', 'D1',
                                   varying=['PD', 'T1'],
                                   fixed=['B1'],
                                   extra={'algo': traits.String(desc="Choose algorithm (l/w/n)", argstr="--algo=%s"),
                                          'iterations': traits.Int(desc='Max iterations for WLLS/NLLS (default 15)', argstr='--its=%d')})

HIFI, HIFISim = base.Command('HIFI', 'qi despot1hifi', 'HIFI',
                             varying=['PD', 'T1', 'B1'],
                             files=['spgr', 'mprage'],
                             extra={'clamp_T1': traits.Float(desc='Clamp T1 between 0 and value', argstr='--clamp=%f')})

DESPOT2, DESPOT2Sim = base.Command('DESPOT2', 'qi despot2', 'D2',
                                   varying=['PD', 'T2'],
                                   fixed=['T1', 'B1'],
                                   extra={'algo': traits.Enum("LLS", "WLS", "NLS", desc="Choose algorithm", argstr="--algo=%d"),
                                          'ellipse': traits.Bool(desc="Data is ellipse geometric solution", argstr='--gs'),
                                          'iterations': traits.Int(desc='Max iterations for WLLS/NLLS (default 15)', argstr='--its=%d'),
                                          'clamp_PD': traits.Float(desc='Clamp PD between 0 and value', argstr='-f %f'),
                                          'clamp_T2': traits.Float(desc='Clamp T2 between 0 and value', argstr='--clampT1=%f')})

FM, FMSim = base.Command('FM', 'qi despot2fm', 'FM',
                         varying=['PD', 'T2', 'f0'],
                         fixed=['B1', 'T1'],
                         extra={'asym': traits.Bool(desc="Fit asymmetric (+/-) off-resonance frequency", argstr='--asym'),
                                'algo': traits.Enum("LLS", "WLS", "NLS", desc="Choose algorithm", argstr="--algo=%d")})

JSR, JSRSim = base.Command('JSR', 'qi jsr', 'JSR',
                           varying=['PD', 'T1', 'T2', 'df0'],
                           fixed=['B1'],
                           files=['spgr', 'ssfp'],
                           extra={'npsi': traits.Int(desc='Number of psi/off-resonance starts', argstr='--npsi=%d')})

Multiecho, MultiechoSim = base.Command('Multiecho',
                                       'qi multiecho',
                                       'ME',
                                       varying=['PD', 'T2'],
                                       extra={'algo': traits.String(desc="Choose algorithm (l/a/n)", argstr="--algo=%s"),
                                              'iterations': traits.Int(desc='Max iterations for WLLS/NLLS (default 15)', argstr='--its=%d'),
                                              'thresh_PD': traits.Float(desc='Only output maps when PD exceeds threshold value', argstr='-t=%f'),
                                              'clamp_T2': traits.Float(desc='Clamp T2 between 0 and value', argstr='-p=%f')})

MPMR2s, MPMR2sSim = base.Command('MPMR2s',
                                 'qi mpm_r2s',
                                 'MPM',
                                 varying=['R2s', 'S0_PDw', 'S0_T1w', 'S0_MTw'],
                                 files=['PDw', 'T1w', 'MTw'])

Ellipse, EllipseSim = base.Command('Ellipse',
                                   'qi ssfp_ellipse',
                                   'ES',
                                   varying=['G', 'a', 'b',
                                            'theta_0', 'phi_rf'],
                                   extra={'algo': traits.String(desc='Choose algorithm (h/d)', argstr='--algo=%s')})

PLANET, PLANETSim = base.Command('PLANET',
                                 'qi planet',
                                 'PLANET',
                                 varying=['PD', 'T1', 'T2'],
                                 fixed=['B1'],
                                 files=['G', 'a', 'b'])


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
