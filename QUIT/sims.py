import requests
import gzip
import nibabel as nib
import numpy as np
from os.path import isfile


def init_brainweb():
    """
    Download the Brainweb Phantoms if they do not already exist in the current directory
    """
    if not isfile('classes.mnc'):
        print('Downloading classes')
        params = {'download_for_real': '[Start Download!]',
                  'do_download_alias': 'phantom_1.0mm_normal_crisp',
                  'format_value': 'minc',
                  'who_name': 'Tobias Wood',
                  'who_institution': 'KCL',
                  'who_email': 'tobias.wood@kcl.ac.uk'}
        response = requests.get(
            url='http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1', params=params)
        minc_file = open('classes.mnc', 'wb')
        minc_file.write(response.content)
    if not isfile('rf20_C.mnc'):
        print('Downloading B1')
        params = {'download_for_real': '[Start Download!]',
                  'do_download_alias': 'rf20_C.mnc.gz',
                  'format_value': 'minc',
                  'zip_value': 'none',
                  'who_name': 'Tobias Wood',
                  'who_institution': 'KCL',
                  'who_email': 'tobias.wood@kcl.ac.uk'}
        response = requests.get(
            url='https://brainweb.bic.mni.mcgill.ca/cgi/brainweb1', params=params)
        b1_file = open('rf20_C.mnc', 'wb')
        b1_file.write(response.content)


def make_phantom(parameters, subsamp=1, B1=True, f0=False):
    # Full class list:
    # 0=Background
    # 1=CSF
    # 2=Grey Matter
    # 3=White Matter
    # 4=Fat
    # 5=Muscle/Skin
    # 6=Skin
    # 7=Skull
    # 8=Glial Matter
    # 9=Connective
    classes = nib.load('classes.mnc')
    class_data = classes.get_data().astype('int32')
    for key, vals in parameters.items():
        data = np.choose(
            class_data[::subsamp, ::subsamp, ::subsamp], np.array(vals)).astype('float32')
        img = nib.nifti1.Nifti1Image(data, affine=classes.affine)
        nib.save(img, key + '.nii.gz')

    if B1:
        b1_img = nib.load('rf20_C.mnc')
        b1_data = b1_img.get_data().astype(
            'float32')[::subsamp, ::subsamp, ::subsamp]
        img = nib.nifti1.Nifti1Image(b1_data, affine=b1_img.affine)
        nib.save(img, 'B1.nii.gz')
    if f0:
        shape = [int(np.ceil(x / subsamp)) for x in class_data.shape]
        f0data = np.tile(
            np.linspace(-60, 60, shape[2]), [shape[0], shape[1], 1])
        img = nib.nifti1.Nifti1Image(f0data, affine=classes.affine)
        nib.save(img, 'f0.nii.gz')
