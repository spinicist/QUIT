import numpy as np
from nipype import Workflow, Node, MapNode, IdentityInterface
from nipype.interfaces.fsl import maths, BET, MCFLIRT, FLIRT, ApplyXFM, ImageMaths, Smooth
from qipype.fitting import LtzWaterMT, LtzWaterMTSim, LtzAmNOE
from qipype.commands import ZSpec
from qipype.utils import PCA, Select


def prep(zfrqs, dummies=0, pca_retain=0, name='CEST_prep'):
    inputnode = Node(IdentityInterface(fields=['zspec_file', 'ref_file']),
                     name='inputnode')
    outputnode = Node(IdentityInterface(fields=['zspec_file', 'f0_map',
                                                'mask_file', 'ref_file',
                                                'DS', 'MT']),
                      name='outputnode')

    prep = Workflow(name=name)

    moco = Node(MCFLIRT(cost='mutualinfo', mean_vol=True),
                name='moco')
    mask = Node(BET(mask=True, no_output=True), name='mask')
    prep.connect([(inputnode, moco, [('zspec_file', 'in_file'), ('ref_file', 'ref_file')]),
                  (moco, mask, [('mean_img', 'in_file')])])

    zspec_norm = Node(ImageMaths(op_string='-div', out_file='zspec.nii.gz'),
                      name='zspec_norm', iterfield=['in_file'])
    if (dummies > 0):
        ref_index = dummies - 1
        zspec_select = Node(Select(volumes=list(range(dummies, len(zfrqs))), out_file='zspec.nii.gz'),
                            name='zspec_select')
        zfrqs = np.array(zfrqs[dummies:])
        prep.connect([(moco, zspec_select, [('out_file', 'in_file')]),
                      (zspec_select, zspec_norm, [('out_file', 'in_file')])])
    else:
        ref_index = 0
        zfrqs = np.array(zfrqs)
        prep.connect([(moco, zspec_norm, [('out_file', 'in_file')])])
    zspec_ref = Node(Select(volumes=[ref_index, ], out_file='reference.nii.gz'),
                     name='zspec_ref')
    prep.connect([(moco, zspec_ref, [('out_file', 'in_file')]),
                  (zspec_ref, zspec_norm, [('out_file', 'in_file2')])])

    f0_indices = (np.abs(zfrqs) > 7) | (np.abs(zfrqs) < 1.1)
    sat_frqs = zfrqs[f0_indices]
    sat_angles = np.repeat(180.0, len(f0_indices))
    f0_select = Node(Select(volumes=np.where(f0_indices)[0].tolist(), out_file='background_zspec.nii.gz'),
                     name='f0_select')
    prep.connect([(zspec_norm, f0_select, [('out_file', 'in_file')])])

    sequence = {'MTSat': {'pulse': {'p1': 0.4,
                                    'p2': 0.3,
                                    'bandwidth': 0.39},
                          'Trf': 0.02,
                          'TR': 4,
                          'FA': 5,
                          'sat_f0': sat_frqs.tolist(),
                          'sat_angle': sat_angles.tolist()}}
    f0_fit = Node(LtzWaterMT(sequence=sequence, verbose=True),
                  name='f0_fit')
    prep.connect([(f0_select, f0_fit, [('out_file', 'in_file')]),
                  (mask, f0_fit, [('mask_file', 'mask_file')])])

    out_frqs = np.sort(zfrqs)
    f0_correct = Node(ZSpec(in_freqs=zfrqs.tolist(), out_freqs=out_frqs.tolist(),
                            verbose=True), name='f0_correct')
    prep.connect([(zspec_norm, f0_correct, [('out_file', 'in_file')]),
                  (f0_fit, f0_correct, [('DS_f0_map', 'f0_map')]),
                  (mask, f0_correct, [('mask_file', 'mask_file')])])

    if pca_retain > 0:
        f0_pca = Node(
            PCA(retain=pca_retain, projections_file='proj.nii.gz'), name='f0_pca')
        prep.connect([(f0_correct, f0_pca, [('out_file', 'in_file')]),
                      (f0_pca, outputnode, [('out_file', 'zspec_file')])
                      ])
    else:
        prep.connect([(f0_correct, outputnode, [('out_file', 'zspec_file')])])

    prep.connect([(moco, outputnode, [('mean_img', 'ref_file')]),
                  (mask, outputnode, [('out_file', 'mask_file')]),
                  (f0_fit, outputnode, [('DS_f0_map', 'f0_map'),
                                        ('DS_A_map', 'DS'),
                                        ('MT_A_map', 'MT')])])

    return (prep, out_frqs)


def cert(zfrqs):
    inputnode = Node(IdentityInterface(fields=['cert_180', 'cert_360', 'mask_file']),
                     name='inputnode')
    outputnode = Node(IdentityInterface(fields=['cert_spectrum', 'cert_amide']),
                      name='outputnode')

    cert_sub = Node(ImageMaths(op_string='-sub', out_file='cert.nii.gz'),
                    name='cert_subtract', iterfield=['in_file'])
    amide_index = (np.abs(zfrqs - 3.5)).argmin()
    amide = Node(Select(volumes=[amide_index], out_file='amide.nii.gz'),
                 name='select_amide')

    cert = Workflow(name='CERT')
    cert.connect([(inputnode, cert_sub, [('cert_360', 'in_file'), ('cert_180', 'in_file2')]),
                  (cert_sub, amide, [('out_file', 'in_file')]),
                  (cert_sub, outputnode, [('out_file', 'cert_spectrum')]),
                  (amide, outputnode, [('out_file', 'cert_amide')])
                  ])

    return cert


def amide_noe(zfrqs, name='Amide_NOE'):
    inputnode = Node(IdentityInterface(fields=['zspec_file', 'mask_file']),
                     name='inputnode')
    outputnode = Node(IdentityInterface(fields=['zspec', 'DS_A', 'DS_fwhm', 'DS_f0', 'MT_A', 'MT_fwhm', 'MT_f0', 'Amide', 'NOE']),
                      name='outputnode')
    # Repeat 2-pool fit
    f0_indices = (np.abs(zfrqs) > 9.9) | (np.abs(zfrqs) < 1.1)
    sequence = {'MTSat': {'pulse': {'p1': 0.4,
                                    'p2': 0.3,
                                    'bandwidth': 0.39},
                          'Trf': 0.02,
                          'TR': 4,
                          'FA': 5,
                          'sat_f0': zfrqs[f0_indices].tolist(),
                          'sat_angle': np.repeat(180.0, len(f0_indices)).tolist()}}
    backg_select = Node(Select(volumes=np.where(f0_indices)[0].tolist(), out_file='bg_zspec.nii.gz'),
                        name='backg_select')
    backg_fit = Node(LtzWaterMT(sequence=sequence, verbose=True),
                     name='backg_fit')
    # Simulate data for all frequencies
    sequence['MTSat']['sat_f0'] = zfrqs.tolist()
    sequence['MTSat']['sat_angle'] = np.repeat(180.0, len(zfrqs)).tolist()
    backg_sim = Node(LtzWaterMTSim(sequence=sequence,
                                   noise=0,
                                   out_file='backg_sim.nii.gz',
                                   verbose=True),
                     name='backg_sim')
    backg_sub = Node(ImageMaths(op_string='-sub',
                                out_file='no_backg_sub.nii.gz'),
                     name='backg_sub', iterfield=['in_file'])

    f0_indices = (np.abs(zfrqs) > 0.99) & (np.abs(zfrqs) < 10.1)
    sequence['MTSat']['sat_f0'] = zfrqs[f0_indices].tolist()
    sequence['MTSat']['sat_angle'] = np.repeat(
        180.0, len(f0_indices)).tolist()
    an_select = Node(Select(volumes=np.where(f0_indices)[0].tolist(), out_file='fg_zspec.nii.gz'),
                     name='an_select')
    an_fit=Node(LtzAmNOE(sequence=sequence,
                           Zref=0.0,
                           additive=True,
                           verbose=True),
                  name='an_pool')

    wf=Workflow(name=name)
    wf.connect([(inputnode, backg_select, [('zspec_file', 'in_file')]),
                (backg_select, backg_fit, [('out_file', 'in_file')]),
                (inputnode, backg_fit, [('mask_file', 'mask_file')]),
                (backg_fit, backg_sim, [('DS_f0_map', 'DS_f0_map'),
                                        ('DS_fwhm_map', 'DS_fwhm_map'),
                                        ('DS_A_map', 'DS_A_map'),
                                        ('MT_f0_map', 'MT_f0_map'),
                                        ('MT_fwhm_map', 'MT_fwhm_map'),
                                        ('MT_A_map', 'MT_A_map')]),
                (backg_sim, backg_sub, [('out_file', 'in_file')]),
                (inputnode, backg_sub, [('zspec_file', 'in_file2')]),
                (backg_sub, an_select, [('out_file', 'in_file')]),
                (an_select, an_fit, [('out_file', 'in_file')]),
                (inputnode, an_fit, [('mask_file', 'mask_file')]),
                (backg_fit, outputnode, [('DS_A_map', 'DS_A'),
                                         ('DS_fwhm_map', 'DS_fwhm'),
                                         ('DS_f0_map', 'DS_f0'),
                                         ('MT_A_map', 'MT_A'),
                                         ('MT_fwhm_map', 'MT_fwhm'),
                                         ('MT_f0_map', 'MT_f0')]),
                (an_fit, outputnode, [('Amide_A_map', 'Amide'),
                                      ('NOE_A_map', 'NOE')])
                ])
    return wf
