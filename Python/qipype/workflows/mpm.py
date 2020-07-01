import numpy as np
from nipype import Workflow, Node, MapNode, IdentityInterface
from nipype.interfaces.fsl import maths, BET, MCFLIRT, ExtractROI, FLIRT, Merge, ApplyXFM, ImageMaths, BinaryMaths, ConvertXFM
import nipype.interfaces.utility as util
from qipype.interfaces.relax import MPMR2s
from qipype.interfaces.mt import MTSat
from qipype.interfaces.b1 import B1Minus
from qipype.interfaces.fsl import ApplyXfm4D


def init_mpm_wf(me_params, mtsat_params):
    inputnode = Node(IdentityInterface(fields=['pdw_file', 't1w_file', 'mtw_file',
                                               'pdw_cal', 't1w_cal', 'mtw_cal',
                                               'b1_map']),
                     name='inputnode')
    outputnode = Node(IdentityInterface(fields=['pd_map', 'r1_map', 'r2s_map', 'mtsat_map']),
                      name='outputnode')

    wf = Workflow(name='Multi-Parametric-Mapping')

    bet = Node(BET(mask=True, no_output=True), name='brain_mask')
    wf.connect([(inputnode, bet, [('t1w_file', 'in_file')])])

    pd_b1minus = Node(B1Minus(), name='pd_b1minus')
    t1_b1minus = Node(B1Minus(), name='t1_b1minus')
    mt_b1minus = Node(B1Minus(), name='mt_b1minus')
    wf.connect([(inputnode, pd_b1minus, [('pdw_cal', 'in_file')]),
                (inputnode, t1_b1minus, [('t1w_cal', 'in_file')]),
                (inputnode, mt_b1minus, [('mtw_cal', 'in_file')])])

    pd_b1m_hires = Node(FLIRT(apply_xfm=True, uses_qform=True), name='pd_b1m_hires')
    t1_b1m_hires = Node(FLIRT(apply_xfm=True, uses_qform=True), name='t1_b1m_hires')
    mt_b1m_hires = Node(FLIRT(apply_xfm=True, uses_qform=True), name='mt_b1m_hires')
    wf.connect([(pd_b1minus, pd_b1m_hires, [('out_file', 'in_file')]),
                (inputnode, pd_b1m_hires, [('pdw_file', 'reference')]),
                (t1_b1minus, t1_b1m_hires, [('out_file', 'in_file')]),
                (inputnode, t1_b1m_hires, [('t1w_file', 'reference')]),
                (mt_b1minus, mt_b1m_hires, [('out_file', 'in_file')]),
                (inputnode, mt_b1m_hires, [('mtw_file', 'reference')])])
    
    pd_cal = Node(ImageMaths(op_string='-div'), name='pd_cal', iterfield=['in_file'])
    t1_cal = Node(ImageMaths(op_string='-div'), name='t1_cal', iterfield=['in_file'])
    mt_cal = Node(ImageMaths(op_string='-div'), name='mt_cal', iterfield=['in_file'])
    wf.connect([(inputnode, pd_cal, [('pdw_file', 'in_file')]),
                (pd_b1m_hires, pd_cal, [('out_file', 'in_file2')]),
                (inputnode, t1_cal, [('t1w_file', 'in_file')]),
                (t1_b1m_hires, t1_cal, [('out_file', 'in_file2')]),
                (inputnode, mt_cal, [('mtw_file', 'in_file')]),
                (mt_b1m_hires, mt_cal, [('out_file', 'in_file2')])])

    t1_reg = Node(FLIRT(uses_qform=True, cost='mutualinfo'), name='t1_reg')
    mt_reg = Node(FLIRT(uses_qform=True, cost='mutualinfo'), name='mt_reg')
    t1_apply = Node(ApplyXfm4D(single_matrix=True), name='t1_apply')
    mt_apply = Node(ApplyXfm4D(single_matrix=True), name='mt_apply')
    wf.connect([(t1_cal, t1_reg, [('out_file', 'in_file')]),
                (pd_cal, t1_reg, [('out_file', 'reference')]),
                (t1_cal, t1_apply, [('out_file', 'in_file')]),
                (pd_cal, t1_apply, [('out_file', 'ref_vol')]),
                (t1_reg, t1_apply, [('out_matrix_file', 'trans_file')]),
                (mt_cal, mt_reg, [('out_file', 'in_file')]),
                (pd_cal, mt_reg, [('out_file', 'reference')]),
                (mt_cal, mt_apply, [('out_file', 'in_file')]),
                (pd_cal, mt_apply, [('out_file', 'ref_vol')]),
                (mt_reg, mt_apply, [('out_matrix_file', 'trans_file')])])

    mpm = Node(MPMR2s(sequence=me_params, verbose=True), name='MPM_R2s')
    mtsat = Node(MTSat(sequence=mtsat_params, verbose=True), name='MPM_MTSat')

    wf.connect([(pd_cal, mpm, [('out_file', 'pdw_file')]),
                (t1_apply, mpm, [('out_file', 't1w_file')]),
                (mt_apply, mpm, [('out_file', 'mtw_file')]),
                (bet, mpm, [('mask_file', 'mask_file')]),
                (mpm, mtsat, [('s0_pdw', 'pdw_file'),
                              ('s0_t1w', 't1w_file'),
                              ('s0_mtw', 'mtw_file')]),
                (inputnode, mtsat, [('b1_map', 'b1_map')]),
                (bet, mtsat, [('mask_file', 'mask_file')]),
                (mpm, outputnode, [('r2s_map', 'r2s_map')]),
                (mtsat, outputnode, [('s0_map', 'pd_map'),
                                     ('r1_map', 'r1_map'),
                                     ('delta_map', 'mtsat_map')])])
    return wf
