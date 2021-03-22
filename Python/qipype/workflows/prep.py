import numpy as np
from nipype import Workflow, Node, MapNode, IdentityInterface
from nipype.interfaces.fsl import maths, BET, MCFLIRT, ExtractROI, FLIRT, Merge, ApplyXFM, ImageMaths, BinaryMaths, ConvertXFM
from nipype.interfaces.ants import Registration, ApplyTransforms
import nipype.interfaces.utility as util
from qipype.utils import Mask, RFProfile, Complex, Filter, CoilCombine
from qipype.fsl import ApplyXfm4D


def COMPOSER(verbose=False, is_bruker=False):
    inputnode = Node(IdentityInterface(
        fields=['in_file', 'ref_file']), name='input')
    outputnode = Node(IdentityInterface(
        fields=['out_file']), name='output')
    wf = Workflow(name='COMPOSER')

    in_mag = Node(Complex(magnitude_out_file='in_mag.nii.gz',
                          verbose=verbose), name='in_magnitude')
    ref_mag = Node(Complex(magnitude_out_file='ref_mag.nii.gz',
                           verbose=verbose), name='ref_magnitude')
    if is_bruker:
        wf.connect([(inputnode, in_mag, [('in_file', 'realimag')])])
        wf.connect([(inputnode, ref_mag, [('ref_file', 'realimag')])])
    else:
        wf.connect([(inputnode, in_mag, [('in_file', 'complex')])])
        wf.connect([(inputnode, ref_mag, [('ref_file', 'complex')])])

    in_mean = Node(maths.MeanImage(), name='in_mean')
    ref_mean = Node(maths.MeanImage(), name='ref_mean')
    wf.connect([(in_mag, in_mean, [('magnitude_out_file', 'in_file')]),
                (ref_mag, ref_mean, [('magnitude_out_file', 'in_file')])])

    register = Node(Registration(dimension=3,
                                 initial_moving_transform_com=1,
                                 transforms=['Rigid'],
                                 metric=['Mattes'],
                                 metric_weight=[1],
                                 transform_parameters=[(0.1,)],
                                 number_of_iterations=[[1000, 500, 250]],
                                 collapse_output_transforms=False,
                                 initialize_transforms_per_stage=False,
                                 radius_or_number_of_bins=[32],
                                 sampling_strategy=['Regular', None],
                                 sampling_percentage=[0.25, None],
                                 convergence_threshold=[1.e-6],
                                 smoothing_sigmas=[[4, 2, 1]],
                                 shrink_factors=[[8, 4, 2]],
                                 sigma_units=['vox'],
                                 output_warped_image=True,
                                 verbose=True), name='register')
    wf.connect([(in_mean, register, [('out_file', 'moving_image')]),
                (ref_mean, register, [('out_file', 'fixed_image')])])

    if is_bruker:
        resample = Node(ApplyTransforms(
            dimension=3, input_image_type=3), name='resample_reference')
        in_x = Node(Complex(complex_out_file='in_x.nii.gz',
                            verbose=verbose), name='in_x')
        ref_x = Node(Complex(complex_out_file='ref_x.nii.gz',
                             verbose=verbose), name='ref_x')
        cc = Node(CoilCombine(), name='cc')
        wf.connect([(inputnode, resample, [('ref_file', 'input_image')]),
                    (in_mean, resample, [('out_file', 'reference_image')]),
                    (register, resample, [
                     ('reverse_transforms', 'transforms')]),
                    (inputnode, in_x, [('in_file', 'realimag')]),
                    (resample, ref_x, [('output_image', 'realimag')]),
                    (in_x, cc, [('complex_out_file', 'in_file')]),
                    (ref_x, cc, [('complex_out_file', 'composer_file')]),
                    (cc, outputnode, [('out_file', 'out_file')])])
    else:
        raise('Not Yet Supported')

    return wf


def init_b1_mcf(rf_pulse=None, scale=150):
    inputnode = Node(IdentityInterface(fields=['2db1map_file', 'ref_file']),
                     name='inputnode')
    outputnode = Node(IdentityInterface(fields=['b1_plus', 'b1_pulse']),
                      name='outputnode')

    b1_b1 = Node(ExtractROI(t_min=0, t_size=1), name='b1_extract_b1')
    b1_filter = Node(Filter(filter_spec='Gauss,3.0'), name='b1_filter')
    b1_mag = Node(ExtractROI(t_min=1, t_size=1), name='b1_extract_mag')

    b1_reg = Node(FLIRT(out_file='b1mag_reg.nii.gz', out_matrix_file='b1mag_reg.mat'),
                  name='b1_reg')
    b1_invert = Node(ConvertXFM(invert_xfm=True), name='b1_invert')
    b1_apply = Node(FLIRT(apply_xfm=True), name='b1_reg_apply')
    b1_scale = Node(ImageMaths(op_string='-div %f' % scale), name='b1_scale')

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
                (b1_scale, outputnode, [('out_file', 'b1_plus')])])
    if rf_pulse:
        b1_rf = Node(
            RFProfile(rf=rf_pulse, out_file='b1_rf.nii.gz'), name='b1_rf')
        wf.connect([(b1_scale, b1_rf, [('out_file', 'in_file')]),
                    (b1_rf, outputnode, [('out_file', 'b1_pulse')])])
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
