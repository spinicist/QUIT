import numpy as np
from nipype import Workflow, Node, MapNode, IdentityInterface
from nipype.interfaces.fsl import maths, BET, MCFLIRT, FLIRT, ApplyXFM, ImageMaths, Smooth
from QUIT.mt import Lorentzian, LorentzianSim, ZSpec
from QUIT.utils import PCA, Select


def init_cest_prep(zfrqs, dummies=0, pca_retain=0, name='CEST_prep'):
    inputnode = Node(IdentityInterface(fields=['zspec_file', 'ref_file']),
                     name='inputnode')
    outputnode = Node(IdentityInterface(fields=['zspec_file', 'f0_map',
                                                'mask_file', 'ref_file',
                                                'DS', 'MT']),
                      name='outputnode')

    moco = Node(MCFLIRT(cost='mutualinfo', mean_vol=True),
                name='moco')
    mask = Node(BET(mask=True, no_output=True), name='mask')

    if (dummies > 0):
        ref_index = dummies - 1
        zspec_select = Node(Select(volumes=range(dummies, len(zfrqs)),
                         name='zspec_select')
        zfrqs=np.array(zfrqs[dummies:])
    else:
        ref_index=0
        zfrqs=np.array(zfrqs)

    zspec_ref=Node(Select(volumes=(ref_index,)),
                     name='zspec_ref')
    zspec_norm=Node(ImageMaths(op_string='-div', out_file='zspec.nii.gz'),
                      name='zspec_norm')

    f0_indices=(np.abs(zfrqs) > 7) | (np.abs(zfrqs) < 1.1)
    sat_frqs=zfrqs[f0_indices]
    sat_angles=np.repeat(180.0, len(f0_indices))
    f0_select=Node(Select(volumes=np.where(f0_indices)[0].tolist()),
                     name='f0_select')
    sequence={'MTSat': {'pulse': {'p1': 0.4,
                                    'p2': 0.3,
                                    'bandwidth': 0.39},
                          'Trf': 0.02,
                          'TR': 4,
                          'FA': 5,
                          'sat_f0': sat_frqs.tolist(),
                          'sat_angle': sat_angles.tolist()}}
    two_pools=[{'name': 'DS',
                  'df0': [0, -2.5, 2.5],
                  'fwhm': [1.0, 1.e-6, 3.0],
                  'A': [0.2, 1.e-3, 1.0],
                  'use_bandwidth': True},
                 {'name': 'MT',
                  'df0': [-2.5, -5.0, -0.5],
                  'fwhm': [50.0, 35.0, 200.0],
                  'A': [0.3, 1.e-3, 1.0]}]
    f0_fit=Node(Lorentzian(sequence=sequence, pools=two_pools, verbose=True),
                  name='f0_fit')

    out_frqs=np.sort(zfrqs)
    f0_correct=Node(ZSpec(in_freqs=zfrqs.tolist(), out_freqs=out_frqs.tolist(),
                            verbose=True), name='f0_correct')

    prep=Workflow(name=name)
    prep.connect([(inputnode, moco, [('zspec_file', 'in_file'), ('ref_file', 'ref_file')]),
                  (moco, zspec_ref, [('out_file', 'in_file')]),
                  (moco, mask, [('mean_img', 'in_file')]),
                  (zspec_ref, zspec_norm, [('roi_file', 'in_file2')]),
                  (zspec_norm, f0_select, [('out_file', 'in_file')]),
                  (f0_select, f0_fit, [('out_file', 'in_file')]),
                  (mask, f0_fit, [('mask_file', 'mask_file')]),
                  (zspec_norm, f0_correct, [('out_file', 'in_file')]),
                  (f0_fit, f0_correct, [('DS_f0', 'f0_map')]),
                  (mask, f0_correct, [('mask_file', 'mask_file')]),
                  (moco, outputnode, [('mean_img', 'ref_file')]),
                  (mask, outputnode, [('out_file', 'mask_file')]),
                  (f0_fit, outputnode, [('DS_f0', 'f0_map'),
                                        ('DS_A', 'DS'),
                                        ('MT_A', 'MT')])
                  ])
    if (dummies > 0):
        prep.connect([(moco, zspec_roi, [('out_file', 'in_file')]),
                      (zspec_roi, zspec_norm, [('roi_file', 'in_file')])])
    else:
        prep.connect([(moco, zspec_norm, [('out_file', 'in_file')])])

    if pca_retain > 0:
        f0_pca=Node(
            PCA(retain=pca_retain, projections_file='proj.nii.gz'), name='f0_pca')
        prep.connect([(f0_correct, f0_pca, [('out_file', 'in_file')]),
                      (f0_pca, outputnode, [('out_file', 'zspec_file')])
                      ])
    else:
        prep.connect([(f0_correct, outputnode, [('out_file', 'zspec_file')])])

    return (prep, out_frqs)


def init_cert(zfrqs):
    inputnode=Node(IdentityInterface(fields=['cert_180', 'cert_360', 'mask_file']),
                     name='inputnode')
    outputnode=Node(IdentityInterface(fields=['cert_spectrum', 'cert_amide']),
                      name='outputnode')

    cert_sub=Node(ImageMaths(op_string='-sub', out_file='cert.nii.gz'),
                    name='cert_subtract')
    amide_index=(np.abs(zfrqs - 3.5)).argmin()
    amide=Node(ExtractROI(t_min=amide_index, t_size=1),
                 name='amide_extract')

    cert=Workflow(name='CERT')
    cert.connect([(inputnode, cert_sub, [('cert_360', 'in_file'), ('cert_180', 'in_file2')]),
                  (cert_sub, amide, [('out_file', 'in_file')]),
                  (cert_sub, outputnode, [('out_file', 'cert_spectrum')]),
                  (amide, outputnode, [('roi_file', 'cert_amide')])
                  ])

    return cert


def init_amide_noe(zfrqs):
    inputnode=Node(IdentityInterface(fields=['zspec_file', 'mask_file']),
                     name='inputnode')
    outputnode=Node(IdentityInterface(fields=['diff_file', 'DS', 'MT', 'Amide', 'NOE']),
                      name='outputnode')
    # Repeat 2-pool fit
    f0_indices=(np.abs(zfrqs) > 9.9) | (np.abs(zfrqs) < 1.1)
    sequence={'MTSat': {'pulse': {'p1': 0.4,
                                    'p2': 0.3,
                                    'bandwidth': 0.39},
                          'Trf': 0.02,
                          'TR': 4,
                          'FA': 5,
                          'sat_f0': zfrqs[f0_indices].tolist(),
                          'sat_angle': np.repeat(180.0, len(f0_indices)).tolist()}}
    two_pools=[{'name': 'DS',
                  'df0': [0, -2.5, 2.5],
                  'fwhm': [1.0, 1.e-6, 3.0],
                  'A': [0.2, 1.e-3, 1.0],
                  'use_bandwidth': True},
                 {'name': 'MT',
                  'df0': [-2.5, -5.0, -0.5],
                  'fwhm': [50.0, 35.0, 200.0],
                  'A': [0.3, 1.e-3, 1.0]}]
    backg_select=Node(Select(volumes=np.where(f0_indices)[0].tolist()),
                        name='backg_select')
    backg_fit=Node(Lorentzian(sequence=sequence,
                                pools=two_pools,
                                verbose=True),
                     name='backg_fit')
    # Simulate data for all frequencies
    sequence['MTSat']['sat_f0']=zfrqs.tolist()
    sequence['MTSat']['sat_angle']=np.repeat(180.0, len(zfrqs)).tolist()
    backg_sim=Node(LorentzianSim(sequence=sequence,
                                   pools=two_pools,
                                   noise=0,
                                   in_file='backg_sim.nii.gz',
                                   verbose=True),
                     name='backg_sim')
    backg_sub=Node(ImageMaths(op_string='-sub',
                                out_file='no_backg_sub.nii.gz'),
                     name='backg_sub')

    an_pools=[{'name': 'Amide',
                 'df0': [3.5, 2.0, 6.0],
                 'fwhm': [2.0, 0.4, 4.0],
                 'A': [0.2, 1.e-3, 0.2],
                 'use_bandwidth': True},
                {'name': 'NOE',
                 'df0': [-4.0, -6.0, -2.0],
                 'fwhm': [2.0, 0.4, 4.0],
                 'A': [0.2, 1.e-3, 0.2],
                 'use_bandwidth': True}]
    f0_indices=(np.abs(zfrqs) > 0.99) & (np.abs(zfrqs) < 10.1)
    sequence['MTSat']['sat_f0']=zfrqs[f0_indices].tolist()
    sequence['MTSat']['sat_angle']=np.repeat(
        180.0, len(f0_indices)).tolist()
    an_select=Node(SelectVols(volumes=np.where(f0_indices)[0].tolist()),
                     name='an_select')

    an_fit=Node(Lorentzian(sequence=sequence,
                             pools=an_pools,
                             Zref=0.0,
                             additive=True,
                             verbose=True),
                  name='an_pool')

    outputnode=Node(IdentityInterface(fields=['zspec',
                                                'f0_map',
                                                'mask',
                                                'DS',
                                                'MT',
                                                'Amide',
                                                'NOE']),
                      name='outputnode')
    wf=Workflow(name='Amide_NOE')
    wf.connect([(inputnode, backg_select, [('zspec_file', 'in_file')]),
                (backg_select, backg_fit, [('out_file', 'in_file')]),
                (inputnode, backg_fit, [('mask_file', 'mask_file')]),
                (backg_fit, backg_sim, [('DS_f0', 'DS_f0'),
                                        ('DS_fwhm', 'DS_fwhm'),
                                        ('DS_A', 'DS_A'),
                                        ('MT_f0', 'MT_f0'),
                                        ('MT_fwhm', 'MT_fwhm'),
                                        ('MT_A', 'MT_A')]),
                (inputnode, backg_sub, [('zspec_file', 'in_file2')]),
                (backg_sim, backg_sub, [('out_file', 'in_file')]),
                (backg_sub, an_select, [('out_file', 'in_file')]),
                (an_select, an_fit, [('out_file', 'in_file')]),
                (inputnode, an_fit, [('mask_file', 'mask_file')]),
                (backg_sub, outputnode, [('out_file', 'diff_file')]),
                (backg_fit, outputnode, [('DS_A', 'DS'), ('MT_A', 'MT')]),
                (an_fit, outputnode, [('Amide_A', 'Amide'),
                                      ('NOE_A', 'NOE')])
                ])
    return wf
