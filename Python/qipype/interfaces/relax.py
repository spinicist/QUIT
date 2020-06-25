#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT relaxometry tools
Requires that the QUIT tools are in your your system path
"""

from os import path, getcwd
from nipype.interfaces.base import TraitedSpec, File, traits
from .. import base


############################ qidespot1 ############################

class DESPOT1InputSpec(base.FitInputSpec):
    # Additional Options
    B1_map = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    algo = traits.String(desc="Choose algorithm (l/w/n)", argstr="--algo=%s")
    iterations = traits.Int(
        desc='Max iterations for WLLS/NLLS (default 15)', argstr='--its=%d')


class DESPOT1(base.FitCommand):
    """
    Run DESPOT1 analysis with qidespot1

    Example
    -------
    >>> from qipype.interfaces.relax import DESPOT1
    >>> seq = {'SPGR': {'TR':5E-3, 'FA':[5,10]} }
    >>> d1 = DESPOT1(sequence=seq, in_file='SPGR.nii.gz')
    >>> d1_res = d1.run()
    >>> print(d1_res.outputs)

    """

    _cmd = 'qi despot1'
    input_spec = DESPOT1InputSpec
    output_spec = base.FitOutputSpec('D1', ['PD', 'T1'])

############################ qidespot1sim ############################


class DESPOT1Sim(base.SimCommand):
    """
    Run DESPOT1 simulation with qidespot1

    Example with parameter dictionary
    -------
    >>> from qipype.interfaces.relax import DESPOT1Sim
    >>> seq = {'SPGR': {'TR':5E-3, 'FA':[5,10]}}
    >>> d1sim = DESPOT1Sim(sequence=seq, PD_map='PD.nii.gz', T1_map='T1.nii.gz', out_file='SPGR.nii.gz')
    >>> d1sim_res = d1.run()
    >>> print(d1sim_res.outputs)
    """

    _cmd = 'qi despot1'
    input_spec = base.SimInputSpec('D1', ['PD', 'T1'], ['B1'])
    output_spec = base.SimOutputSpec('D1')

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
        desc='Clamp T1 between 0 and value', argstr='--clamp=%f')


class HIFI(base.FitCommand):
    """
    Calculate T1 & B1 map with the DESPOT1-HIFI method
    """

    _cmd = 'qi despot1hifi'
    input_spec = HIFIInputSpec
    output_spec = base.FitOutputSpec('HIFI', ['PD', 'T1', 'B1'])

############################ qidespot1hifisim ############################


class HIFISim(base.SimCommand):
    """
    Simulate SPGR/FLASH and MPRAGE images using DESPOT1-HIFI model
    """

    _cmd = 'qi despot1hifi'
    sim_files = ['spgr', 'mprage']
    input_spec = base.SimInputSpec(
        'HIFI', ['PD', 'T1', 'B1'], out_files=sim_files)
    output_spec = base.SimOutputSpec('HIFI', sim_files)

############################ qidespot2 ############################


class DESPOT2InputSpec(base.FitInputSpec):
    # Inputs
    T1_map = File(exists=True, argstr='--T1=%s', mandatory=True,
                  position=-2, desc='Path to T1 map')

    # Options
    B1_map = File(desc='B1 map (ratio) file', argstr='--B1=%s')
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


class DESPOT2(base.FitCommand):
    """
    Run DESPOT2 analysis with qidespot2
    """

    _cmd = 'qi despot2'
    input_spec = DESPOT2InputSpec
    output_spec = base.FitOutputSpec('D2', ['PD', 'T2'])

############################ qidespot2sim ############################


class DESPOT2Sim(base.SimCommand):
    """
    Run DESPOT2 simulation
    """

    _cmd = 'qi despot2'
    input_spec = base.SimInputSpec('D2',
                                   varying=['PD', 'T2'],
                                   fixed=['B1', 'T1'],
                                   extras={'ellipse': traits.Bool(
                                           desc="Data is ellipse geometric solution", argstr='--gs')})
    output_spec = base.SimOutputSpec('D2')

############################ qidespot2fm ############################


class FMInputSpec(base.FitInputSpec):
    # Inputs
    T1_map = File(exists=True, argstr='--T1=%s',
                  mandatory=True, desc='Path to T1 map')

    # Options
    B1_map = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    asym = traits.Bool(
        desc="Fit asymmetric (+/-) off-resonance frequency", argstr='--asym')
    algo = traits.Enum("LLS", "WLS", "NLS",
                       desc="Choose algorithm", argstr="--algo=%d")


class FM(base.FitCommand):
    """
    Run DESPOT2-FM analysis
    """

    _cmd = 'qi despot2fm'
    input_spec = FMInputSpec
    output_spec = base.FitOutputSpec('FM', ['PD', 'T2', 'f0'])


class FMSim(base.SimCommand):
    """
    Run DESPOT2-FM simulation
    """

    _cmd = 'qi despot2fm'
    _param_files = ['PD', 'T2', 'f0']
    input_spec = base.SimInputSpec('FM',
                                   varying=['PD', 'T2', 'f0'],
                                   fixed=['B1', 'T1'],
                                   extras={'asym': traits.Bool(
                                       desc="Fit asymmetric (+/-) off-resonance frequency", argstr='--asym')})
    output_spec = base.SimOutputSpec('FM')

############################ qi_jsr ############################


class JSRInputSpec(base.InputSpec):
    # Inputs
    spgr_file = File(exists=True, argstr='%s', mandatory=True,
                     position=-2, desc='Path to SPGR data')

    ssfp_file = File(exists=True, argstr='%s', mandatory=True,
                     position=-1, desc='Path to SSFP data')

    # Options
    B1_map = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    residuals = traits.Bool(
        desc='Write out residuals for each data-point', argstr='--resids')
    npsi = traits.Int(
        desc='Number of psi/off-resonance starts', argstr='--npsi=%d')


class JSR(base.FitCommand):
    """
    Calculate T1 &T2 map with Joint System Relaxometry

    """

    _cmd = 'qi jsr'
    input_spec = JSRInputSpec
    output_spec = base.FitOutputSpec('JSR', ['PD', 'T1', 'T2', 'df0'])

############################ qimcdespot ############################
# Status: Everything is there but not tested. Very much unsupported.


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


class QIMCDespot(base.FitCommand):
    """
    Interace for qimcdespot

    Example 1
    -------
    >>> from qipype.interfaces.relax import QiMcDespot
    >>> interface = QiMcDespot(prefix='nipype_', param_file='mcdespot_params.json')
    """

    _cmd = 'qi mcdespot'
    input_spec = QIMCDespotInputSpec
    output_spec = base.FitOutputSpec(
        '3C', ['T1_m', 'T2_m', 'T1_ie', 'T2_ie', 'T1_csf', 'T2_csf', 'tau_m', 'f_m', 'f_csf', 'f0', 'B1'])

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


class Multiecho(base.FitCommand):
    """
    Run a multi-echo fit
    """

    _cmd = 'qi multiecho'
    input_spec = MultiechoInputSpec
    output_spec = base.FitOutputSpec('ME', ['PD', 'T2'])


class MultiechoSim(base.SimCommand):
    """
    Run a multi-echo simulation
    """

    _cmd = 'qi multiecho'
    input_spec = base.SimInputSpec('ME', ['PD', 'T2'])
    output_spec = base.SimOutputSpec('ME')

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


class MPMR2s(base.FitCommand):
    """
    Runs qi_mpm_r2s
    """

    _cmd = 'qi mpm_r2s'
    input_spec = MPMR2sInputSpec
    output_spec = base.FitOutputSpec(
        'MPM', ['R2s', 'S0_PDw', 'S0_T1w', 'S0_MTw'])


############################ qi_ssfp_elipse ############################


class EllipseInputSpec(base.FitInputSpec):
    # Additional Options
    algo = traits.String(desc='Choose algorithm (h/d)', argstr='--algo=%s')


class Ellipse(base.FitCommand):
    """
    Fit an ellipse to SSFP data
    """

    _cmd = 'qi ssfp_ellipse'
    input_spec = EllipseInputSpec
    output_spec = base.FitOutputSpec(
        'ES', ['G', 'a', 'b', 'theta_0', 'phi_rf'])


class EllipseSim(base.SimCommand):
    """
    Simulate SSFP data from ellipse parameters
    """

    _cmd = 'qi ssfp_ellipse'
    input_spec = base.SimInputSpec('ES', ['G', 'a', 'b', 'theta_0', 'phi_rf'])
    output_spec = base.SimOutputSpec('ES')

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
    B1_map = File(desc='B1 map (ratio) file', argstr='--B1=%s')
    residuals = traits.Bool(
        desc='Write out residuals for each data-point', argstr='--resids')


class PLANET(base.FitCommand):
    """
    Calculate T1/T2 from ellipse parameters

    """

    _cmd = 'qi planet'
    input_spec = PLANETInputSpec
    output_spec = base.FitOutputSpec('PLANET', ['PD', 'T1', 'T2'])


class PLANETSim(base.SimCommand):
    """
    Simulate ellipse parameters from T1/T2

    """

    _cmd = 'qi planet'
    input_spec = base.SimInputSpec('PLANET',
                                   varying=['PD', 'T1', 'T2'],
                                   fixed=['B1'],
                                   out_files=['G', 'a', 'b'])
    output_spec = base.SimOutputSpec('PLANET', ['G', 'a', 'b'])

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['G_file'] = path.abspath(self.inputs.G_file)
        outputs['a_file'] = path.abspath(self.inputs.a_file)
        outputs['b_file'] = path.abspath(self.inputs.b_file)
        return outputs
