import numpy as np
from nipype import Workflow, Node, MapNode, IdentityInterface
from nipype.interfaces.fsl import maths, BET, MCFLIRT, ExtractROI, FLIRT, Merge, ApplyXFM, ImageMaths, BinaryMaths, ConvertXFM
import nipype.interfaces.utility as util
from qipype.fitting import MPMR2s, MTSat
from qipype.commands import B1Minus
from qipype.fsl import ApplyXfm4D
from qipype.utils import Select


def init_mpm_b1_wf(me_params, mtsat_params):
    inputnode = Node(IdentityInterface(fields=['PDw_file', 'T1w_file', 'MTw_file',
                                               'PDw_cal', 'T1w_cal', 'MTw_cal',
                                               'B1_map']),
                     name='inputnode')
    outputnode = Node(IdentityInterface(fields=['PD_map', 'R1_map', 'R2s_map', 'mtsat_map']),
                      name='outputnode')

    wf = Workflow(name='Multi-Parametric-Mapping')

    bet = Node(BET(mask=True, no_output=True), name='brain_mask')
    wf.connect([(inputnode, bet, [('T1w_file', 'in_file')])])

    PD_B1minus = Node(B1Minus(), name='PD_B1minus')
    t1_B1minus = Node(B1Minus(), name='t1_B1minus')
    mt_B1minus = Node(B1Minus(), name='mt_B1minus')
    wf.connect([(inputnode, PD_B1minus, [('PDw_cal', 'in_file')]),
                (inputnode, t1_B1minus, [('T1w_cal', 'in_file')]),
                (inputnode, mt_B1minus, [('MTw_cal', 'in_file')])])

    pd0 = Node(Select(volumes=[0, ], out_file='pd0.nii.gz'), name='pd0')
    t10 = Node(Select(volumes=[0, ], out_file='t10.nii.gz'), name='t10')
    mt0 = Node(Select(volumes=[0, ], out_file='mt0.nii.gz'), name='mt0')
    wf.connect([(inputnode, pd0, [('PDw_file', 'in_file')]),
                (inputnode, t10, [('T1w_file', 'in_file')]),
                (inputnode, mt0, [('MTw_file', 'in_file')])])

    PD_B1m_hires = Node(
        FLIRT(apply_xfm=True, uses_qform=True), name='PD_B1m_hires')
    t1_B1m_hires = Node(
        FLIRT(apply_xfm=True, uses_qform=True), name='t1_B1m_hires')
    mt_B1m_hires = Node(
        FLIRT(apply_xfm=True, uses_qform=True), name='mt_B1m_hires')
    wf.connect([(PD_B1minus, PD_B1m_hires, [('out_file', 'in_file')]),
                (pd0, PD_B1m_hires, [('out_file', 'reference')]),
                (t1_B1minus, t1_B1m_hires, [('out_file', 'in_file')]),
                (t10, t1_B1m_hires, [('out_file', 'reference')]),
                (mt_B1minus, mt_B1m_hires, [('out_file', 'in_file')]),
                (mt0, mt_B1m_hires, [('out_file', 'reference')])])

    pd_cal = Node(ImageMaths(op_string='-div'),
                  name='pd_cal', iterfield=['in_file'])
    t1_cal = Node(ImageMaths(op_string='-div'),
                  name='t1_cal', iterfield=['in_file'])
    mt_cal = Node(ImageMaths(op_string='-div'),
                  name='mt_cal', iterfield=['in_file'])
    wf.connect([(inputnode, pd_cal, [('PDw_file', 'in_file')]),
                (PD_B1m_hires, pd_cal, [('out_file', 'in_file2')]),
                (inputnode, t1_cal, [('T1w_file', 'in_file')]),
                (t1_B1m_hires, t1_cal, [('out_file', 'in_file2')]),
                (inputnode, mt_cal, [('MTw_file', 'in_file')]),
                (mt_B1m_hires, mt_cal, [('out_file', 'in_file2')])])

    t1_reg = Node(FLIRT(uses_qform=True, cost='mutualinfo'), name='t1_reg')
    mt_reg = Node(FLIRT(uses_qform=True, cost='mutualinfo'), name='mt_reg')
    t1_apply = Node(ApplyXfm4D(single_matrix=True), name='t1_apply')
    mt_apply = Node(ApplyXfm4D(single_matrix=True), name='mt_apply')
    wf.connect([(t10, t1_reg, [('out_file', 'in_file')]),
                (pd0, t1_reg, [('out_file', 'reference')]),
                (t1_cal, t1_apply, [('out_file', 'in_file')]),
                (pd0, t1_apply, [('out_file', 'ref_vol')]),
                (t1_reg, t1_apply, [('out_matrix_file', 'trans_file')]),
                (mt0, mt_reg, [('out_file', 'in_file')]),
                (pd0, mt_reg, [('out_file', 'reference')]),
                (mt_cal, mt_apply, [('out_file', 'in_file')]),
                (pd0, mt_apply, [('out_file', 'ref_vol')]),
                (mt_reg, mt_apply, [('out_matrix_file', 'trans_file')])])

    mpm = Node(MPMR2s(sequence=me_params, verbose=True), name='MPM_R2s')
    mtsat = Node(MTSat(sequence=mtsat_params, verbose=True), name='MPM_MTSat')

    wf.connect([(pd_cal, mpm, [('out_file', 'PDw_file')]),
                (t1_apply, mpm, [('out_file', 'T1w_file')]),
                (mt_apply, mpm, [('out_file', 'MTw_file')]),
                (bet, mpm, [('mask_file', 'mask_file')]),
                (mpm, mtsat, [('S0_PDw_map', 'PDw_file'),
                              ('S0_T1w_map', 'T1w_file'),
                              ('S0_MTw_map', 'MTw_file')]),
                (inputnode, mtsat, [('B1_map', 'B1_map')]),
                (bet, mtsat, [('mask_file', 'mask_file')]),
                (mpm, outputnode, [('R2s_map', 'R2s_map')]),
                (mtsat, outputnode, [('PD_map', 'PD_map'),
                                     ('R1_map', 'R1_map'),
                                     ('delta_map', 'mtsat_map')])])
    return wf


def init_mpm_wf(me_params, mtsat_params):
    inputnode = Node(IdentityInterface(fields=['PDw_file', 'T1w_file', 'MTw_file']),
                     name='inputnode')
    outputnode = Node(IdentityInterface(fields=['PD_map', 'R1_map', 'R2s_map', 'mtsat_map']),
                      name='outputnode')

    wf = Workflow(name='Multi-Parametric-Mapping')

    bet = Node(BET(mask=True, no_output=True), name='brain_mask')
    wf.connect([(inputnode, bet, [('T1w_file', 'in_file')])])

    pd0 = Node(Select(volumes=[0, ], out_file='pd0.nii.gz'), name='pd0')
    t10 = Node(Select(volumes=[0, ], out_file='t10.nii.gz'), name='t10')
    mt0 = Node(Select(volumes=[0, ], out_file='mt0.nii.gz'), name='mt0')
    wf.connect([(inputnode, pd0, [('PDw_file', 'in_file')]),
                (inputnode, t10, [('T1w_file', 'in_file')]),
                (inputnode, mt0, [('MTw_file', 'in_file')])])

    t1_reg = Node(FLIRT(uses_qform=True, cost='mutualinfo'), name='t1_reg')
    mt_reg = Node(FLIRT(uses_qform=True, cost='mutualinfo'), name='mt_reg')
    t1_apply = Node(ApplyXfm4D(single_matrix=True), name='t1_apply')
    mt_apply = Node(ApplyXfm4D(single_matrix=True), name='mt_apply')
    wf.connect([(t10, t1_reg, [('out_file', 'in_file')]),
                (pd0, t1_reg, [('out_file', 'reference')]),
                (inputnode, t1_apply, [('T1w_file', 'in_file')]),
                (pd0, t1_apply, [('out_file', 'ref_vol')]),
                (t1_reg, t1_apply, [('out_matrix_file', 'trans_file')]),
                (mt0, mt_reg, [('out_file', 'in_file')]),
                (pd0, mt_reg, [('out_file', 'reference')]),
                (inputnode, mt_apply, [('MTw_file', 'in_file')]),
                (pd0, mt_apply, [('out_file', 'ref_vol')]),
                (mt_reg, mt_apply, [('out_matrix_file', 'trans_file')])])

    mpm = Node(MPMR2s(sequence=me_params, verbose=True), name='MPM_R2s')
    mtsat = Node(MTSat(sequence=mtsat_params, verbose=True), name='MPM_MTSat')

    wf.connect([(inputnode, mpm, [('PDw_file', 'PDw_file')]),
                (t1_apply, mpm, [('out_file', 'T1w_file')]),
                (mt_apply, mpm, [('out_file', 'MTw_file')]),
                (bet, mpm, [('mask_file', 'mask_file')]),
                (mpm, mtsat, [('S0_PDw_map', 'PDw_file'),
                              ('S0_T1w_map', 'T1w_file'),
                              ('S0_MTw_map', 'MTw_file')]),
                (inputnode, mtsat, [('B1_map', 'B1_map')]),
                (bet, mtsat, [('mask_file', 'mask_file')]),
                (mpm, outputnode, [('R2s_map', 'R2s_map')]),
                (mtsat, outputnode, [('PD_map', 'PD_map'),
                                     ('R1_map', 'R1_map'),
                                     ('delta_map', 'mtsat_map')])])
    return wf
