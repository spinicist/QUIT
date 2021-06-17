from pathlib import Path
from os import chdir
import unittest
import numpy as np
from nipype.interfaces.base import CommandLine
from qipype.fitting import LtzWaterMT, LtzWaterMTSim, qMT, qMTSim
from qipype.commands import NewImage, Diff, ZSpec, Lineshape

vb = True
CommandLine.terminal_output = 'allatonce'


class MT(unittest.TestCase):
    def setUp(self):
        Path('testdata').mkdir(exist_ok=True)
        chdir('testdata')

    def tearDown(self):
        chdir('../')

    def test_ltz(self):
        sat_f0 = [*np.linspace(-40, 40,
                               12).squeeze().tolist(), -0.5, -0.25, 0, 0.25, 0.5]
        sequence = {'MTSat': {'pulse': {'p1': 0.4,
                                        'p2': 0.3,
                                        'bandwidth': 0.39},
                              'TR': 4,
                              'Trf': 0.02,
                              'FA': 5,
                              'sat_f0': sat_f0,
                              'sat_angle': np.repeat(180.0, 17).squeeze().tolist()}}
        lorentz_file = 'lorentz_sim.nii.gz'
        img_sz = [32, 32, 32]
        noise = 0.001

        NewImage(out_file='PD.nii.gz', verbose=vb, img_size=img_sz,
                 fill=1.0).run()
        NewImage(out_file='f0.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=0, grad_vals=(-0.25, 0.25)).run()

        NewImage(out_file='fwhm.nii.gz', verbose=vb, img_size=img_sz,
                 fill=1.8).run()
        NewImage(out_file='A.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=2, grad_vals=(0.3, 0.4)).run()

        NewImage(out_file='MTf0.nii.gz', verbose=vb, img_size=img_sz,
                 fill=-2.2).run()
        NewImage(out_file='MTfwhm.nii.gz', verbose=vb, img_size=img_sz,
                 fill=90).run()
        NewImage(out_file='MTA.nii.gz', verbose=vb, img_size=img_sz,
                 fill=0.4).run()

        LtzWaterMTSim(sequence=sequence, out_file=lorentz_file,
                      noise=noise, verbose=vb,
                      DS_f0_map='f0.nii.gz',
                      DS_fwhm_map='fwhm.nii.gz',
                      DS_A_map='A.nii.gz',
                      MT_fwhm_map='MTfwhm.nii.gz',
                      MT_f0_map='f0.nii.gz',
                      MT_A_map='MTA.nii.gz').run()
        LtzWaterMT(sequence=sequence, in_file=lorentz_file, verbose=vb).run()

        diff_fwhm = Diff(in_file='LTZ_DS_fwhm.nii.gz', baseline='fwhm.nii.gz',
                         noise=noise, verbose=vb).run()
        diff_A = Diff(in_file='LTZ_DS_A.nii.gz', baseline='A.nii.gz',
                      noise=noise, verbose=vb).run()

        diff_MTfwhm = Diff(in_file='LTZ_MT_fwhm.nii.gz', baseline='MTfwhm.nii.gz',
                           noise=noise, verbose=vb).run()
        diff_MTA = Diff(in_file='LTZ_MT_A.nii.gz', baseline='MTA.nii.gz',
                        noise=noise, verbose=vb).run()

        self.assertLessEqual(diff_fwhm.outputs.out_diff, 15)
        self.assertLessEqual(diff_A.outputs.out_diff, 300)
        self.assertLessEqual(diff_MTfwhm.outputs.out_diff, 15)
        self.assertLessEqual(diff_MTA.outputs.out_diff, 300)

    def test_qMT(self):
        qmt = {'MTSat': {
            'TR': 0.032,
            'Trf': 0.020,
            'FA': 5,
            'sat_f0': [1000, 1000, 2236, 2236, 5000, 5000, 11180, 11180, 250000, 250000],
            'sat_angle': [750, 360, 750, 360, 750, 360, 750, 360, 750, 360],
            'pulse': {'name': 'Gauss', 'p1': 0.416, 'p2': 0.295, 'bandwidth': 200}
        }
        }
        qmt_file = 'qmt_sim.nii.gz'
        t1app = 'qmt_t1app.nii.gz'
        img_sz = [32, 32, 1]
        noise = 0.001

        lineshape_file = '_qmt_lineshape.json'

        Lineshape(out_file=lineshape_file, lineshape='SuperLorentzian',
                  frq_start=500, frq_space=500, frq_count=150).run()

        NewImage(out_file='M0_f.nii.gz', verbose=vb, img_size=img_sz,
                 fill=1.0).run()
        NewImage(out_file='f_b.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=1, grad_vals=(0.01, 0.15)).run()
        NewImage(out_file='T2_b.nii.gz', verbose=vb, img_size=img_sz,
                 fill=12e-6).run()
        NewImage(out_file='T2_f.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=0, grad_vals=(0.05, 0.15)).run()
        NewImage(out_file='k.nii.gz', verbose=vb, img_size=img_sz,
                 fill=4.0).run()
        NewImage(out_file=t1app, verbose=vb, img_size=img_sz,
                 fill=1.0).run()

        qMTSim(sequence=qmt, out_file=qmt_file, T1_map=t1app, lineshape=lineshape_file,
               noise=noise, verbose=vb,
               M0_f_map='M0_f.nii.gz',
               f_b_map='f_b.nii.gz',
               T2_b_map='T2_b.nii.gz',
               T2_f_map='T2_f.nii.gz',
               k_map='k.nii.gz').run()
        qMT(sequence=qmt, in_file=qmt_file, T1_map=t1app,
            lineshape=lineshape_file, verbose=vb).run()

        diff_PD = Diff(in_file='QMT_M0_f.nii.gz', baseline='M0_f.nii.gz',
                       noise=noise, verbose=vb).run()
        diff_f_b = Diff(in_file='QMT_f_b.nii.gz', baseline='f_b.nii.gz',
                        noise=noise, verbose=vb).run()
        diff_k = Diff(in_file='QMT_k.nii.gz', baseline='k.nii.gz',
                      noise=noise, verbose=vb).run()
        self.assertLessEqual(diff_PD.outputs.out_diff, 10)
        self.assertLessEqual(diff_f_b.outputs.out_diff, 30)
        self.assertLessEqual(diff_k.outputs.out_diff, 35)

    def test_ZSpec(self):
        NewImage(out_file='zspec_linear.nii.gz', verbose=vb, img_size=[8, 8, 8, 4],
                 grad_dim=3, grad_vals=(-3, 3)).run()
        NewImage(out_file='zspec_zero.nii.gz', verbose=vb,
                 img_size=[8, 8, 8], fill=0).run()
        ZSpec(in_file='zspec_linear.nii.gz',
              in_freqs=[-3, -1, 1, 3],
              out_freqs=[0],
              verbose=vb).run()
        diff_zero = Diff(in_file='zspec_linear_interp.nii.gz', abs_diff=True,
                         baseline='zspec_zero.nii.gz', noise=1, verbose=vb).run()
        self.assertLessEqual(diff_zero.outputs.out_diff, 0.01)


if __name__ == '__main__':
    unittest.main()
