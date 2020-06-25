import numpy as np
from nipype import Workflow, Node, MapNode, IdentityInterface
from nipype.interfaces.fsl import maths, BET, MCFLIRT, ExtractROI, FLIRT, Merge, ApplyXFM, ImageMaths, BinaryMaths, ConvertXFM
import nipype.interfaces.utility as util
from quit.interfaces.relax import MPMR2s
from quit.interfaces.mt import MTSat
from .interfaces import ApplyXfm4D


def init_mpm_wf(me_params, mtsat_params):
    inputnode = Node(IdentityInterface(fields=['pdw_file', 't1w_file', 'mtw_file',
                                               'pdw_cal', 't1w_cal', 'mtw_cal']),
                     name='inputnode')
    outputnode = Node(IdentityInterface(fields=['pd_map', 'r1_map', 'r2s_map', 'mtsat_map']),
                      name='outputnode')

    bet = Node(BET(mask=True, no_output=True), name='brain_mask')
    mpm = Node(MPMR2s(sequence=me_params, verbose=True), name='MPM_R2s')
    mtsat = Node(MTSat(sequence=mtsat_params, verbose=True), name='MPM_MTSat')

    wf = Workflow(name='Multi-Parametric-Mapping')
    wf.connect([(inputnode, bet, [('t1w_file', 'in_file')]),
                (inputnode, mpm, [('pdw_file', 'pdw_file'),
                                  ('t1w_file', 't1w_file'),
                                  ('mtw_file', 'mtw_file')]),
                (bet, mpm, [('mask_file', 'mask_file')]),
                (mpm, mtsat, [('s0_pdw', 'pdw_file'),
                              ('s0_t1w', 't1w_file'),
                              ('s0_mtw', 'mtw_file')]),
                (bet, mtsat, [('mask_file', 'mask_file')]),
                (mpm, outputnode, [('r2s_map', 'r2s_map')]),
                (mtsat, outputnode, [('s0_map', 'pd_map'),
                                     ('r1_map', 'r1_map'),
                                     ('delta_map', 'mtsat_map')])])
    return wf
