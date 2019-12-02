import numpy as np
from nipype import Workflow, Node, MapNode, IdentityInterface
from nipype.interfaces.fsl import maths, BET, MCFLIRT, ExtractROI, FLIRT, Merge, ApplyXFM, ImageMaths, BinaryMaths, ConvertXFM
import nipype.interfaces.utility as util
from QUIT.interfaces.utils import Mask, RFProfile, Complex, Filter
from .interfaces import ApplyXfm4D


def init_b1_mcf(rf_pulse, scale=150):
    inputnode = Node(IdentityInterface(fields=['2db1map_file', 'ref_file']),
                     name='inputnode')
    outputnode = Node(IdentityInterface(fields=['b1_map']),
                      name='outputnode')

    b1_b1 = Node(ExtractROI(t_min=0, t_size=1), name='b1_extract_b1')
    b1_filter = Node(Filter(filter_spec='Gauss,3.0'), name='b1_filter')
    b1_mag = Node(ExtractROI(t_min=1, t_size=1), name='b1_extract_mag')

    b1_reg = Node(FLIRT(out_file='b1mag_reg.nii.gz', out_matrix_file='b1mag_reg.mat'),
                  name='b1_reg')
    b1_invert = Node(ConvertXFM(invert_xfm=True), name='b1_invert')
    b1_apply = Node(FLIRT(apply_xfm=True), name='b1_reg_apply')
    b1_scale = Node(ImageMaths(op_string='-div %f' % scale), name='b1_scale')
    b1_rf = Node(RFProfile(rf=rf_pulse, out_file='b1_rf.nii.gz'), name='b1_rf')

    wf = Workflow(name='b1_prep')
    wf.connect([(inputnode, b1_b1, [('2db1map_file', 'in_file')]),
                (inputnode, b1_mag, [('2db1map_file', 'in_file')]),
                (inputnode, b1_reg, [('ref_file', 'in_file')]),
                (inputnode, b1_apply, [('ref_file', 'reference')]),
                (b1_mag, b1_reg, [('roi_file', 'reference')]),
                (b1_reg, b1_invert, [('out_matrix_file', 'in_file')]),
                (b1_invert, b1_apply, [('out_file', 'in_matrix_file')]),
                (b1_b1, b1_filter, [('roi_file', 'in_file')]),
                (b1_filter, b1_apply, [('out_file', 'in_file')]),
                (b1_apply, b1_scale, [('out_file', 'in_file')]),
                (b1_scale, b1_rf, [('out_file', 'in_file')]),
                (b1_rf, outputnode, [('out_file', 'b1_map')])])
    return wf


def init_complex_mcf(name='', ref=False, fix_ge=True, negate=True):
    inputnode = Node(IdentityInterface(fields=['real_file', 'imag_file', 'ref_file']),
                     name='inputnode')
    outputnode = Node(IdentityInterface(fields=['x_file', 'ref_file', 'mask_file']),
                      name='outputnode')

    ri = Node(Complex(fix_ge=fix_ge, negate=negate,
                      magnitude_out_file=name+'_mag.nii.gz',
                      real_out_file=name+'_r.nii.gz',
                      imag_out_file=name+'_i.nii.gz'),
              name='ri_'+name)
    moco = Node(MCFLIRT(mean_vol=not ref, save_mats=True), name='moco_'+name)
    apply_r = Node(ApplyXfm4D(four_digit=True), name='apply_r' + name)
    apply_i = Node(ApplyXfm4D(four_digit=True), name='apply_i' + name)
    x = Node(Complex(complex_out_file=name+'_x.nii.gz'), name='x_'+name)
    f = Node(Filter(complex_in=True, complex_out=True,
                    filter_spec='Tukey'), name='filter_'+name)

    wf = Workflow(name='prep_' + name)
    wf.connect([(inputnode, ri, [('real_file', 'real')]),
                (inputnode, ri, [('imag_file', 'imag')]),
                (ri, moco, [('magnitude_out_file', 'in_file')]),
                (ri, apply_r, [('real_out_file', 'in_file')]),
                (ri, apply_i, [('imag_out_file', 'in_file')]),
                (moco, apply_r, [('mat_dir', 'trans_dir')]),
                (moco, apply_i, [('mat_dir', 'trans_dir')]),
                (apply_r, x, [('out_file', 'real')]),
                (apply_i, x, [('out_file', 'imag')]),
                (x, f, [('complex_out_file', 'in_file')]),
                (f, outputnode, [('out_file', 'x_file')])])
    if not ref:
        mask = Node(BET(mask=True, no_output=True), name='mask')
        wf.connect([(moco, mask, [('mean_img', 'in_file')]),
                    (moco, apply_r, [('mean_img', 'ref_vol')]),
                    (moco, apply_i, [('mean_img', 'ref_vol')]),
                    (moco, outputnode, [('mean_img', 'ref_file')]),
                    (mask, outputnode, [('mask_file', 'mask_file')])])
    else:
        wf.connect([(inputnode, moco, [('ref_file', 'ref_file')]),
                    (inputnode, apply_r, [('ref_file', 'ref_vol')]),
                    (inputnode, apply_i, [('ref_file', 'ref_vol')])])

    return wf
