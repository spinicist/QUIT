import unittest
from QUIT.base import BaseCommand
from QUIT.core import NewImage, Diff
from QUIT.mt import Lorentzian, LorentzianSim, Lineshape, qMT, qMTSim, ZSpec

vb = True
BaseCommand.terminal_output = 'allatonce'


class MT(unittest.TestCase):
    def test_lorentzian(self):
        lorentz = {'Lorentz': {'sat_f0': [-5, -4, -3, -2, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5]},
                   'PDFile': 'PD.nii.gz',
                   'f0File': 'f0.nii.gz',
                   'fwhmFile': 'fwhm.nii.gz',
                   'AFile': 'A.nii.gz'}
        lorentz_file = 'lorentz_sim.nii.gz'
        img_sz = [32, 32, 32]
        noise = 0.001

        NewImage(out_file='PD.nii.gz', verbose=vb, img_size=img_sz,
                 fill=1.0).run()
        NewImage(out_file='f0.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=0, grad_vals=(-0.5, 0.5)).run()
        NewImage(out_file='fwhm.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=1, grad_vals=(0.5, 5)).run()
        NewImage(out_file='A.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=2, grad_vals=(0.1, 1)).run()
        LorentzianSim(param_dict=lorentz, in_file=lorentz_file,
                      noise=noise, verbose=vb).run()
        Lorentzian(param_dict=lorentz, in_file=lorentz_file, verbose=vb).run()

        diff_PD = Diff(in_file='LTZ_PD.nii.gz', baseline='PD.nii.gz',
                       noise=noise, verbose=vb).run()
        diff_f0 = Diff(in_file='LTZ_f0.nii.gz', baseline='f0.nii.gz',
                       noise=noise, verbose=vb).run()
        diff_fwhm = Diff(in_file='LTZ_fwhm.nii.gz', baseline='fwhm.nii.gz',
                         noise=noise, verbose=vb).run()
        diff_A = Diff(in_file='LTZ_A.nii.gz', baseline='A.nii.gz',
                      noise=noise, verbose=vb).run()
        self.assertLessEqual(diff_PD.outputs.out_diff, 0.6)
        self.assertLessEqual(diff_f0.outputs.out_diff, 50)
        self.assertLessEqual(diff_fwhm.outputs.out_diff, 5)
        self.assertLessEqual(diff_A.outputs.out_diff, 2)

    def test_qMT(self):
        qmt = {'MTSat': {'sat_f0': [1000, 3000, 5000, 7000, 9000, 1000, 3000, 5000, 7000, 9000],
                         'sat_angle': [360, 360, 360, 360, 360, 720, 720, 720, 720, 720],
                         'FA': 5,
                         'TR': 0.03,
                         'pulse': {"name": "Gauss", "Trf": 0.015, "p1": 0.416, "p2": 0.295}},
               'PDFile': 'PD.nii.gz',
                         'T1_fFile': 'T1_f.nii.gz',
                         'T2_fFile': 'T2_f.nii.gz',
                         'T2_bFile': 'T2_b.nii.gz',
                         'k_bfFile': 'k_bf.nii.gz',
                         'f_bFile': 'f_b.nii.gz',
                         'f0File': 'f0.nii.gz',
                         'B1File': 'B1.nii.gz'}
        qmt_file = 'qmt_sim.nii.gz'
        t1app = 'qmt_t1app.nii.gz'
        img_sz = [32, 32, 32]
        noise = 0.001

        lineshape_file = 'qmt_lineshape.json'

        Lineshape(out_file=lineshape_file, lineshape='SuperLorentzian',
                  frq_start=500, frq_space=500, frq_count=150).run()

        NewImage(out_file='PD.nii.gz', verbose=vb, img_size=img_sz,
                 fill=1.0).run()
        NewImage(out_file='T1_f.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=0, grad_vals=(0.5, 1.5)).run()
        NewImage(out_file='T2_f.nii.gz', verbose=vb, img_size=img_sz,
                 fill=0.1).run()
        NewImage(out_file='T2_b.nii.gz', verbose=vb, img_size=img_sz,
                 fill=12e-6).run()
        NewImage(out_file='k_bf.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=1, grad_vals=(1.0, 5.0)).run()
        NewImage(out_file='f_b.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=2, grad_vals=(0.05, 0.2)).run()

        qMTSim(param_dict=qmt, in_file=qmt_file, t1_map=t1app, lineshape=lineshape_file,
               noise=noise, verbose=vb).run()
        qMT(param_dict=qmt, in_file=qmt_file, t1_map=t1app,
            lineshape=lineshape_file, verbose=vb).run()

        diff_T1_f = Diff(in_file='QMT_T1_f.nii.gz', baseline='T1_f.nii.gz',
                         noise=noise, verbose=vb).run()
        diff_k_bf = Diff(in_file='QMT_k_bf.nii.gz', baseline='k_bf.nii.gz',
                         noise=noise, verbose=vb).run()
        diff_f_b = Diff(in_file='QMT_f_b.nii.gz', baseline='f_b.nii.gz',
                        noise=noise, verbose=vb).run()
        self.assertLessEqual(diff_T1_f.outputs.out_diff, 32)
        self.assertLessEqual(diff_k_bf.outputs.out_diff, 32)
        self.assertLessEqual(diff_f_b.outputs.out_diff, 32)

    def test_ZSpec(self):
        zspec = {'input_freqs': [-3, -1, 1, 3],
                 'output_freqs': [0]}
        NewImage(out_file='zspec_linear.nii.gz', verbose=vb, img_size=[8, 8, 8, 4],
                 grad_dim=3, grad_vals=(-3, 3)).run()
        NewImage(out_file='zspec_zero.nii.gz', verbose=vb,
                 img_size=[8, 8, 8], fill=0).run()
        ZSpec(in_file='zspec_linear.nii.gz',
              param_dict=zspec, verbose=vb).run()
        diff_zero = Diff(in_file='zspec_linear_interp.nii.gz', abs_diff=True,
                         baseline='zspec_zero.nii.gz', noise=1, verbose=vb).run()
        self.assertLessEqual(diff_zero.outputs.out_diff, 0.01)


if __name__ == '__main__':
    unittest.main()
