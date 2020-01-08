#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT relaxometry utilities
Requires that the QUIT tools are in your your system path
"""

from os import path, getcwd
from nipype.interfaces.base import TraitedSpec, File, traits
from .. import base as QI


############################ MUPA ############################


class MUPAInputSpec(QI.FitInputSpec):
    # Inputs - none

    # Options - none. Yet
    pass


class MUPAOutputSpec(TraitedSpec):
    pd_map = File('MUPA_PD.nii.gz', desc="Path to PD map", usedefault=True)
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
    # Inputs

    # Options
    pass


class MUPASim(QI.SimCommand):
    """
    Run MUPA simulation

    """

    _cmd = 'qi mupa'
    _param_files = ['PD', 'T1', 'T2']
    input_spec = MUPASimInputSpec
    output_spec = QI.SimOutputSpec

############################ MUPA-MT ############################


class MUPAMTInputSpec(QI.FitInputSpec):
    # Inputs - none

    lineshape = File(desc='Path to linshape file', argstr='--lineshape=%s')


class MUPAMTOutputSpec(TraitedSpec):
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
    input_spec = MUPAMTInputSpec
    output_spec = MUPAMTOutputSpec


class MUPAMTSimInputSpec(QI.SimInputSpec):
    # Inputs

    lineshape = File(desc='Path to linshape file', argstr='--lineshape=%s')


class MUPAMTSim(QI.SimCommand):
    """
    Run MUPA simulation

    """

    _cmd = 'qi mupa --mt'
    _param_files = ['M0_f', 'T1_f', 'T2_f', 'M0_b', 'B1']
    input_spec = MUPAMTSimInputSpec
    output_spec = QI.SimOutputSpec


############################ MT ############################


class MTInputSpec(QI.FitInputSpec):
    # Inputs - none

    lineshape = File(desc='Path to linshape file', argstr='--lineshape=%s')


class MTOutputSpec(TraitedSpec):
    M0_f_map = File('MT_M0_f.nii.gz',
                    desc='Path to free pool M0 map', usedefault=True)
    M0_b_map = File('MT_M0_b.nii.gz',
                    desc='Path to bound pool M0 map', usedefault=True)
    t1_f_map = File('MT_T1_f.nii.gz',
                    desc='Path to T1 map', usedefault=True)
    t2_f_map = File('MT_T2_f.nii.gz',
                    desc='Path to T2 map', usedefault=True)
    t2_b_map = File('MT_T2_b.nii.gz',
                    desc='Path to bound pool T2 map', usedefault=True)
    k_map = File('MT_k.nii.gz', desc='Path to exchange rate map',
                 usedefault=True)
    f_b_map = File('MT_f_b.nii.gz',
                   desc='Path to f_b map', usedefault=True)
    b1_map = File('MT_B1.nii.gz',
                  desc='Path to B1 map', usedefault=True)
    rmse_map = File('MT_rmse.nii.gz',
                    desc="Path to residual map", usedefault=True)


class MT(QI.FitCommand):
    """
    Run  MT Analysis

    """

    _cmd = 'qi rufis-mt'
    input_spec = MTInputSpec
    output_spec = MTOutputSpec


class MTSimInputSpec(QI.SimInputSpec):
    # Inputs

    lineshape = File(desc='Path to linshape file', argstr='--lineshape=%s')


class MTSim(QI.SimCommand):
    """
    Run  simulation

    """

    _cmd = 'qi rufis-mt'
    _param_files = ['M0_f', 'T1_f', 'T2_f', 'M0_b', 'T2_b', 'k', 'B1']
    input_spec = MTSimInputSpec
    output_spec = QI.SimOutputSpec
