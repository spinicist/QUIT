#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT relaxometry utilities
Requires that the QUIT tools are in your your system path
"""

from os import path, getcwd
from nipype.interfaces.base import TraitedSpec, DynamicTraitedSpec, File, traits
from .. import base as QI


############################ MUPA ############################


class MUPAInputSpec(QI.FitInputSpec):
    # Inputs - none

    # Options - none. Yet
    pass


class MUPAOutputSpec(TraitedSpec):
    pd_map = File('MUPA_M0.nii.gz', desc="Path to M0 map", usedefault=True)
    t1_map = File('MUPA_T1.nii.gz', desc="Path to T1 map", usedefault=True)
    t2_map = File('MUPA_T2.nii.gz', desc="Path to T2 map", usedefault=True)
    rmse_map = File('MUPA_rmse.nii.gz',
                    desc="Path to residual map", usedefault=True)


class MUPA(QI.FitCommand):
    """
    Run MUPA Analysis

    """

    _cmd = 'qi mupa'
    input_spec = MUPAInputSpec
    output_spec = MUPAOutputSpec


class MUPASimInputSpec(QI.SimInputSpec):
    # Options
    pass


class MUPASim(QI.SimCommand):
    """
    Run MUPA simulation

    """

    _cmd = 'qi mupa'
    _param_files = ['M0', 'T1', 'T2']
    input_spec = MUPASimInputSpec
    output_spec = QI.SimOutputSpec

############################ MUPA-B1 ############################


class MUPAB1InputSpec(QI.FitInputSpec):
    # Inputs - none

    # Options - none. Yet
    pass


class MUPAB1OutputSpec(TraitedSpec):
    pd_map = File('MUPAB1_M0.nii.gz', desc="Path to M0 map", usedefault=True)
    t1_map = File('MUPAB1_T1.nii.gz', desc="Path to T1 map", usedefault=True)
    t2_map = File('MUPAB1_T2.nii.gz', desc="Path to T2 map", usedefault=True)
    b1_map = File('MUPAB1_B1.nii.gz', desc="Path to B1 map", usedefault=True)
    rmse_map = File('MUPAB1_rmse.nii.gz',
                    desc="Path to residual map", usedefault=True)


class MUPAB1(QI.FitCommand):
    """
    Run MUPA Analysis

    """

    _cmd = 'qi mupa --B1'
    input_spec = MUPAInputSpec
    output_spec = MUPAOutputSpec


class MUPAB1SimInputSpec(QI.SimInputSpec):
    # Options
    pass


class MUPAB1Sim(QI.SimCommand):
    """
    Run MUPA simulation

    """

    _cmd = 'qi mupa'
    _param_files = ['M0', 'T1', 'T2', 'B1']
    input_spec = MUPASimInputSpec
    output_spec = QI.SimOutputSpec

############################ MUPA-MT ############################


class MUPASteadyStateInputSpec(QI.FitInputSpec):
    # Inputs - none
    pass


class MUPASteadyStateOutputSpec(TraitedSpec):
    M0_f_map = File('MUPAMT_M0_f.nii.gz',
                    desc='Path to free pool M0 map', usedefault=True)
    M0_b_map = File('MUPAMT_M0_b.nii.gz',
                    desc='Path to bound pool M0 map', usedefault=True)
    t1_f_map = File('MUPAMT_T1_f.nii.gz',
                    desc='Path to T1 map', usedefault=True)
    t2_f_map = File('MUPAMT_T2_f.nii.gz',
                    desc='Path to T2 map', usedefault=True)
    f_b_map = File('MUPAMT_f_b.nii.gz',
                   desc='Path to f_b map', usedefault=True)
    b1_map = File('MUPAMT_B1.nii.gz',
                  desc='Path to B1 map', usedefault=True)
    rmse_map = File('MUPA_rmse.nii.gz',
                    desc="Path to residual map", usedefault=True)


class MUPAMT(QI.FitCommand):
    """
    Run MUPA MT Analysis

    """

    _cmd = 'qi mupa --mt'
    input_spec = MUPASteadyStateInputSpec
    output_spec = MUPASteadyStateOutputSpec


class MUPAMTSimInputSpec(QI.SimInputSpec):
    # Inputs
    pass


class MUPAMTSim(QI.SimCommand):
    """
    Run MUPA simulation

    """

    _cmd = 'qi mupa --mt'
    _param_files = ['M0_f', 'M0_b', 'T1_f', 'T2_f', 'B1']
    input_spec = MUPAMTSimInputSpec
    output_spec = QI.SimOutputSpec


############################ Steady-State ############################


class SteadyStateInputSpec(QI.FitInputSpec):
    fitT2 = traits.Bool(desc='Fit T2 model', argstr='--T2')


class SteadyStateOutputSpec(DynamicTraitedSpec):
    pass


class SteadyState(QI.FitCommand):
    """
    Run  Steady-State Analysis

    """

    _cmd = 'qi rufis-ss'
    _prefix = 'SS'
    input_spec = SteadyStateInputSpec
    output_spec = SteadyStateOutputSpec

    def __init__(self, **kwargs):
        if 'fitT2' in kwargs and kwargs['fitT2']:
            self._param_files = ['M0', 'T1', 'T2', 'f0', 'B1']
        else:
            self._param_files = ['M0', 'T1', 'B1']
        super(SteadyState, self).__init__(**kwargs)


class SteadyStateSimInputSpec(QI.SimInputSpec):
    fitT2 = traits.Bool(desc='Fit T2 model', argstr='--T2')


class SteadyStateSim(QI.SimCommand):
    """
    Run  simulation

    """

    _cmd = 'qi rufis-ss'
    input_spec = SteadyStateSimInputSpec
    output_spec = QI.SimOutputSpec

    def __init__(self, **kwargs):
        if 'fitT2' in kwargs and kwargs['fitT2']:
            self._param_files = ['M0', 'T1', 'T2', 'f0', 'B1']
        else:
            self._param_files = ['M0', 'T1', 'B1']
        super(SteadyStateSim, self).__init__(**kwargs)
