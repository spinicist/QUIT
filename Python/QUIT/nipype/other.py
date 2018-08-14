#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Implementation of special nipype interfaces for QUIT
"""

from nipype.interfaces.base import BaseInterfaceInputSpec, BaseInterface, File, TraitedSpec, traits
import nibabel as nib
import numpy as np
import os
from nipype.interfaces import utility as util

def merge_images(img_flist, w_flist, out_img_fname):
    out_img = np.zeros_like(nib.load(img_flist[0]).get_data())
    # Make weighted sum of images
    for i in range(len(img_flist)):
        I = nib.load(img_flist[i]).get_data()
        w = nib.load(w_flist[i]).get_data()
        for t in range(np.size(I,3)):
            out_img[:,:,:,t] = out_img[:,:,:,t] + I[:,:,:,t]*w

    # Save output img
    nii = nib.Nifti1Image(out_img, nib.load(img_flist[0]).affine)
    nib.save(nii, out_img_fname)

class MergeImagesInputSpec(BaseInterfaceInputSpec):
    wm_img = File(exists=True, mandatory=True, desc='WM image')
    gm_img = File(exists=True, mandatory=True, desc='WM image')
    csf_img = File(exists=True, mandatory=True, desc='WM image')
    wm_weight = File(exists=True, mandatory=True, desc='WM weight')
    gm_weight = File(exists=True, mandatory=True, desc='GM weight')
    csf_weight = File(exists=True, mandatory=True, desc='CSF weight')
    out_file = File(mandatory=True, desc='the output image') # Do not set exists=True !!

class MergeImagesOutputSpec(TraitedSpec):
    out_file = File(desc='the output image')

class MergeImages(BaseInterface):
    input_spec = MergeImagesInputSpec
    output_spec = MergeImagesOutputSpec

    def _run_interface(self, runtime):
        # Call our python code here:
        imgs = [self.inputs.wm_img, self.inputs.gm_img, self.inputs.csf_img]
        ws = [self.inputs.wm_weight, self.inputs.gm_weight, self.inputs.csf_weight]
        merge_images(imgs, ws, self.inputs.out_file)
        # And we are done
        return runtime

    def _list_outputs(self):
        return {'out_file': os.path.abspath(self.inputs.out_file)}

###################################################################