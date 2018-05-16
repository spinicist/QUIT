#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of a simulation framework for QUIT. Including downloading data
"""

from tqdm import tqdm
import requests
import math
import os
import zipfile
import shutil
import nibabel as nib
import numpy as np
import json

class Phantom(object):
    """
    Base class for phantom object
    """

    def __init__(self, name, path):
        self.path = path
        self.name = name
        self.sim_files = {}
        self.all_files = {}

class BrainPhantom(Phantom):
    """
    Class for quantitative 3D brain phantom
    """

    def __init__(self, name, path):
        super().__init__(name, path)
        self.tissue_params = {
            'WM':{'T1':None, 'T2':None, 'T2s':None, 'PD':None},
            'WM2':None,
            'WM3':None,
            'GM':{'T1':None, 'T2':None, 'T2s':None, 'PD':None},
            'CSF':{'T1':None, 'T2':None, 'T2s':None, 'PD':None}
        }
        self.brain = None
        self.affine = None
    
    def make_brain(self):
        """
        Create a standard phantom object.
        
        Run the standard routine for creating the quantitative T1/T2/PD maps from
        the ICBM data and the object tissue parameters

        We want a dictionary, and corresponding files, with the following links

        - Tissue class (WM/GM/CSF):
            - T1
            - T2
            - T2s
            - PD
        
        Filenames as: WM_T1.nii.gz (etc).

        Run one simulation for each tissue class and add the resulting images
        """
    
        # Get directory where the module is placed
        package_directory = os.path.dirname(os.path.abspath(__file__))
        p = os.path.join(package_directory, 'sim_data')

        masks = {'CSF':'csf.nii', 'GM':'gm.nii', 'WM':'wm.nii'}
        mask_img = {}
        for k in masks.keys():
            mask_img[k] = nib.load(os.path.join(p, masks[k])).get_data()
        
        # Get nii info for affine to save images later
        nii = nib.load(os.path.join(p,masks['WM']))
        self.affine = nii.affine

        # Make
        self.brain = {'WM':{}, 'GM':{}, 'CSF':{}}
        for tissue in self.brain.keys():
            for param in ['T1','T2','T2s','PD','w']:
                self.brain[tissue][param] = np.zeros_like(mask_img['WM'])

        # Make the images
        for tissue in ['WM','GM','CSF']:
            for param in ['T1','T2','T2s','PD']:
                bin_mask = np.where(mask_img[tissue]>0,1,0)
                self.brain[tissue][param] = bin_mask * self.tissue_params[tissue][param]
            
            # Make an image that corresponds to the tissue weight
            self.brain[tissue]['w'] = mask_img[tissue]

    def save_phantom(self):
        """
        Save brain as .nii files in specified directory and populate the 
        file path dictionary. Also save tissue params as .json

        return: file path dictionary

        -> Check if folder already exists
        """

        # Make folder for files
        full_base_path = os.path.abspath(self.path)
        path_ok = False
        temp_name = self.name
        i = 1
        while path_ok is not True:
            if os.path.exists(os.path.join(full_base_path, temp_name)):
                temp_name = '{}_{}'.format(self.name,i)
                i += 1
            else:
                path_ok = True
        
        brain_path = os.path.join(full_base_path, temp_name)
        os.mkdir(brain_path)

        for tissue in self.brain.keys():
            self.all_files[tissue] = {}
            self.sim_files[tissue] = {}
            for param in ['T1', 'T2', 'T2s', 'PD', 'w']:
                fname = '{}_{}.nii.gz'.format(tissue, param)
                nii = nib.Nifti1Image(self.brain[tissue][param], self.affine)
                full_fname = os.path.join(brain_path,fname)
                nib.save(nii, full_fname)
                # We don't want the weight for the qisignal simulation
                self.all_files[tissue][param] = full_fname
                if param != 'w':
                    self.sim_files[tissue][param] = full_fname
            
            for param in ['f0', 'B1']:
                self.sim_files[tissue][param] = ""

        # Save json file for each tissue
        for tissue in self.brain.keys():
            json_fname = os.path.join(brain_path, '{}_params.json'.format(tissue))
            with open(json_fname,'w') as f:
                json.dump(self.sim_files[tissue], f)

        return self.all_files

    def load_tissues(self, nc=1, field_strength=3):
        """
        Create a default tissue model for the given field strength and number of 
        tissue pools
            nc: Number of tissue compartments in WM. 1=Standard, 2=IE water and Myelin water, 3=IE water, Myelin water, solid myelin
            field_strength: B0 field strength.
        """
        if nc>1:
            print('Only one tissue class is currently supported')
            return
        
        # Params from https://www.researchgate.net/publication/4371528_2D_and_3D_Shepp-logan_phantom_standards_for_MRI
        self.tissue_params['WM']['T1'] = 0.583*field_strength**(0.382)
        self.tissue_params['WM']['T2'] = 0.08
        self.tissue_params['WM']['T2s'] = 0.08
        self.tissue_params['WM']['PD'] = 0.617

        self.tissue_params['GM']['T1'] = 0.857*field_strength**(0.376)
        self.tissue_params['GM']['T2'] = 0.1
        self.tissue_params['GM']['T2s'] = 0.1
        self.tissue_params['GM']['PD'] = 0.745

        self.tissue_params['CSF']['T1'] = 4.2
        self.tissue_params['CSF']['T2'] = 1.99
        self.tissue_params['CSF']['T2s'] = 1.99
        self.tissue_params['CSF']['PD'] = 1

class SLPhantom(Phantom):
    """
    Class for Quantitative 2D Shep-Logan phantom
    """

    def __init__(self, name, path):
        super().__init__(name, path)


def parse_bw_tissue_params():
    package_directory = os.path.dirname(os.path.abspath(__file__))
    fname = os.path.join(package_directory, 'sim_data', 'tissue_mr_parameters.txt')
    f = open(fname, 'r')
    d = {}
    in_entry = False
    ck = None
    for line in f:
        if '#' not in line:
            if 'Tissue Name' in line:
                in_entry = True
                ck = line.split(':')[-1].lstrip().rstrip()
                d[ck] = {}
            else:
                line = line.split(':')
                param = line[0].split('(')[0].lstrip().rstrip()
                if '*' in param:
                    param = param.replace('*', 's')
                val = float(line[-1].lstrip().rstrip())
                d[ck][param] = val

    dd = {}
    for k in d.keys():
        new_k = int(d[k]['Tissue Label'])
        dd[new_k] = d[k].copy()
        dd[new_k].pop('Tissue Label')
        dd[new_k]['Tissue Name'] = k

    return dd

def download_file(url, fname):
	# Streaming, so we can iterate over the response.
	r = requests.get(url, stream=True)

	# Total size in bytes.
	total_size = int(r.headers.get('content-length', 0)); 
	block_size = 1024
	wrote = 0 
	with open(fname, 'wb') as f:
		for data in tqdm(r.iter_content(block_size), desc='Dowloading %s'%fname, total=math.ceil(total_size//block_size) , unit='KB', unit_scale=True, leave=False):
			wrote = wrote  + len(data)
			f.write(data)
	if total_size != 0 and wrote != total_size:
		print("ERROR, something went wrong")  

def download_simdata():
    down_data = {
        # 'SPM':{
        # 	'url':'https://www.fil.ion.ucl.ac.uk/spm/toolbox/TPM/BlaiottaTPM.nii',
        # 	'fname':'BlaiottaTPM.nii'
        # },
        # 'colin27':{
        # 	'url':'http://packages.bic.mni.mcgill.ca/mni-models/colin27/mni_colin27_2008_fuzzy_nifti.zip',
        # 	'fname':'mni_colin27_2008_fuzzy_nifti.zip'	
        # },
        # 'ICBM2009a':{
        # 	'url':"http://www.bic.mni.mcgill.ca/~vfonov/icbm/2009/mni_icbm152_nlin_sym_09a_nifti.zip",
        # 	'fname':"mni_icbm152_nlin_sym_09a_nifti.zip"
        # },
        'ICBM2009c':{
            'url':"http://www.bic.mni.mcgill.ca/~vfonov/icbm/2009/mni_icbm152_nlin_sym_09c_nifti.zip",
            'fname':"mni_icbm152_nlin_sym_09c_nifti.zip"
        },
        'BrainWeb Tissue Params':{
            'url':'http://brainweb.bic.mni.mcgill.ca/tissue_mr_parameters.txt',
            'fname':'tissue_mr_parameters.txt'
        }
    }

    print("----- Downloading Sample Data -----")

    try:
        os.mkdir('sim_data')
    except FileExistsError:
        print('Folder exist already')

    for dset in down_data.keys():
        print('File information:')
        print("URL: "+down_data[dset]['url'])
        print("File: "+down_data[dset]['fname'])
        fpath = os.path.join('sim_data',down_data[dset]['fname'])
        if not os.path.isfile(fpath):
            download_file(down_data[dset]['url'], fpath)
        else:
            print('File already exist, no download required')
        print('---------------------------------------------------')

    fnames = [
        # 'mni_icbm152_nlin_sym_09a_nifti.zip',
        # 'mni_colin27_2008_fuzzy_nifti.zip',
        'mni_icbm152_nlin_sym_09c_nifti.zip'
        ]

    print('------ Unzipping data ------')
    fnames = [os.path.join('sim_data', x) for x in fnames]
    for f in fnames:
        output_fname = f.split('.')[0]
        if not os.path.isfile(output_fname):
            print ('-> Unzipping %s' % f)
            with zipfile.ZipFile(f,'r') as zip_ref:
                zip_ref.extractall(output_fname)

    move_icbm152_09c_data()

def move_icbm152_09c_data():
    data_path = 'sim_data/mni_icbm152_nlin_sym_09c_nifti/mni_icbm152_nlin_sym_09c'
    dest_path = 'sim_data'
    print('-> Moving files')
    os.rename(os.path.join(data_path,'mni_icbm152_csf_tal_nlin_sym_09c.nii'),os.path.join(dest_path, 'csf.nii'))
    os.rename(os.path.join(data_path,'mni_icbm152_wm_tal_nlin_sym_09c.nii'), os.path.join(dest_path, 'wm.nii'))
    os.rename(os.path.join(data_path,'mni_icbm152_gm_tal_nlin_sym_09c.nii'), os.path.join(dest_path, 'gm.nii'))
    os.rename(os.path.join(data_path,'mni_icbm152_t1_tal_nlin_sym_09c_mask.nii'), os.path.join(dest_path, 'mask.nii'))
    os.rename(os.path.join(data_path,'mni_icbm152_t1_tal_nlin_sym_09c.nii'), os.path.join(dest_path, 't1w.nii'))
    os.rename(os.path.join(data_path,'mni_icbm152_t2_tal_nlin_sym_09c.nii'), os.path.join(dest_path, 't2w.nii'))
    os.rename(os.path.join(data_path,'mni_icbm152_pd_tal_nlin_sym_09c.nii'), os.path.join(dest_path, 'pdw.nii'))
    print('-> Removing temp files')
    shutil.rmtree('sim_data/mni_icbm152_nlin_sym_09c_nifti')
    os.remove('sim_data/mni_icbm152_nlin_sym_09c_nifti.zip')

def make_discrete_model(tissue_params):
    print('--- Make discrete model ---')
    minc = nib.load('phantom_1.0mm_normal_crisp.mnc.gz')
    img = minc.get_data()

    T1_map = np.zeros_like(img)
    T2_map = np.zeros_like(img)
    T2s_map = np.zeros_like(img)
    PD_map = np.zeros_like(img)

    for k in tissue_params.keys():
        tclass = int(k)
        T1_map = T1_map + np.where(img==tclass, tissue_params[k]['T1'], 0)
        T2_map = T2_map + np.where(img==tclass, tissue_params[k]['T2'], 0)
        T2s_map = T2s_map + np.where(img==tclass, tissue_params[k]['T2s'], 0)
        PD_map = PD_map + np.where(img==tclass, tissue_params[k]['PD'], 0)

    print('-> Saving T1 map')
    T1_nii = nib.Nifti1Image(T1_map, minc.affine)
    nib.save(T1_nii, 'T1map_disc.nii.gz')

    print('-> Saving T2 map')
    T2_nii = nib.Nifti1Image(T2_map, minc.affine)
    nib.save(T2_nii, 'T2map_disc.nii.gz')

    print('-> Saving T2s map')
    T2s_nii = nib.Nifti1Image(T2s_map, minc.affine)
    nib.save(T2s_nii, 'T2smap_disc.nii.gz')

    print('-> Saving PD map')
    PD_nii = nib.Nifti1Image(PD_map, minc.affine)
    nib.save(PD_nii, 'PDmap_disc.nii.gz')

def make_fussy_model(tissue_params):
    print('--- Making Fuzzy model ---')
    d = {1:'phantom_1.0mm_normal_csf.mnc.gz',
     2:'phantom_1.0mm_normal_gry.mnc.gz',
     3:'phantom_1.0mm_normal_wht.mnc.gz',
     4:'phantom_1.0mm_normal_fat.mnc.gz',
     5:'phantom_1.0mm_normal_m+s.mnc.gz',
     6:'phantom_1.0mm_normal_skn.mnc.gz',
     7:'phantom_1.0mm_normal_skl.mnc.gz',
     8:'phantom_1.0mm_normal_gli.mnc.gz',
     9:'phantom_1.0mm_normal_mit.mnc.gz'}

    
    bw_img = {}
    for k in d.keys():
        bw_img[k] = nib.load(d[k]).get_data()
    minc = nib.load(d[1])

    T1_map = np.zeros_like(bw_img[1])
    T2_map = np.zeros_like(bw_img[1])
    T2s_map = np.zeros_like(bw_img[1])
    PD_map = np.zeros_like(bw_img[1])

    print('-> Generating quantitative maps')
    # Make a weighted tissue map
    for k in range(1,10):
        T1_map = T1_map + bw_img[k]*tp[str(k)]['T1']
        T2_map = T2_map + bw_img[k]*tp[str(k)]['T2']
        T2s_map = T2s_map + bw_img[k]*tp[str(k)]['T2s']
        PD_map = PD_map + bw_img[k]*tp[str(k)]['PD']

    print('-> Saving T1 map')
    T1_nii = nib.Nifti1Image(T1_map, minc.affine)
    nib.save(T1_nii, 'T1map_fuz.nii.gz')

    print('-> Saving T2 map')
    T2_nii = nib.Nifti1Image(T2_map, minc.affine)
    nib.save(T2_nii, 'T2map_fuz.nii.gz')

    print('-> Saving T2s map')
    T2s_nii = nib.Nifti1Image(T2s_map, minc.affine)
    nib.save(T2s_nii, 'T2smap_fuz.nii.gz')

    print('-> Saving PD map')
    PD_nii = nib.Nifti1Image(PD_map, minc.affine)
    nib.save(PD_nii, 'PDmap_fuz.nii.gz')

def make_icbm_09c_model(tissue_params):
    print('--- Making Fuzzy model ---')
    d = {1:'csf.nii', 2:'gm.nii',3:'wm.nii'}
    
    bw_img = {}
    package_directory = os.path.dirname(os.path.abspath(__file__))
    p = os.path.join(package_directory, 'sim_data')
    for k in d.keys():
        bw_img[k] = nib.load(os.path.join(p,d[k])).get_data()
    
    nii = nib.load(os.path.join(p,d[1]))
    tp = tissue_params
    T1_map = np.zeros_like(bw_img[1])
    T2_map = np.zeros_like(bw_img[1])
    T2s_map = np.zeros_like(bw_img[1])
    PD_map = np.zeros_like(bw_img[1])

    print('-> Generating quantitative maps')
    # Make a weighted tissue map
    for k in range(1,4):
        T1_map = T1_map + bw_img[k]*(tp[k]['T1'])
        T2_map = T2_map + bw_img[k]*(tp[k]['T2'])
        T2s_map = T2s_map + bw_img[k]*(tp[k]['T2s'])
        PD_map = PD_map + bw_img[k]*tp[k]['PD']

    print('-> Saving T1 map')
    T1_nii = nib.Nifti1Image(T1_map, nii.affine)
    nib.save(T1_nii, os.path.join(package_directory, 'sim_data','T1map_icbm_09c.nii.gz'))

    print('-> Saving T2 map')
    T2_nii = nib.Nifti1Image(T2_map, nii.affine)
    nib.save(T2_nii, os.path.join(package_directory, 'sim_data','T2map_icbm_09c.nii.gz'))

    print('-> Saving T2s map')
    T2s_nii = nib.Nifti1Image(T2s_map, nii.affine)
    nib.save(T2s_nii, os.path.join(package_directory, 'sim_data','T2smap_icbm_09c.nii.gz'))

    print('-> Saving PD map')
    PD_nii = nib.Nifti1Image(PD_map, nii.affine)
    nib.save(PD_nii, os.path.join(package_directory, 'sim_data','PDmap_fuz.nii.gz'))
