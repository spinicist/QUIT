from pathlib import Path
from os import chdir
import unittest
from nipype.interfaces.base import CommandLine
from qipype.commands import NewImage, Diff
from qipype.parmesan import SteadyStateT1, SteadyStateT1Sim

vb = True
CommandLine.terminal_output = 'allatonce'


class Parmesan(unittest.TestCase):
    def setUp(self):
        Path('testdata').mkdir(exist_ok=True)
        chdir('testdata')

    def tearDown(self):
        chdir('../')

    def test_steadystate(self):
        seq = {
            "TR": 0.002,
            "Tramp": 0.01,
            "Tspoil": 0.05,
            "FA": [2, 8, 8],
            "Trf": 28e-6,
            "spokes_per_seg": 48,
            "groups_per_seg": [1, 1, 1],
            "prep_FA": [0, 0, 360],
            "prep_df": [0, 0, 0],
            "prep_Trf": 0.002,
            "prep_p1": 1.,
            "prep_p2": 1.
        }
        sim_file = 'sim_ss.nii.gz'
        img_sz = [32, 32, 32]
        noise = 0.001

        NewImage(img_size=img_sz, grad_dim=0, grad_vals=(0.8, 1.0),
                 out_file='M0.nii.gz', verbose=vb).run()
        NewImage(img_size=img_sz, grad_dim=1, grad_vals=(0.8, 1.4),
                 out_file='T1.nii.gz', verbose=vb).run()
        NewImage(img_size=img_sz, grad_dim=2, grad_vals=(0.9, 1.1),
                 out_file='B1.nii.gz', verbose=vb).run()

        SteadyStateT1Sim(sequence=seq, out_file=sim_file,
                         M0_map='M0.nii.gz', T1_map='T1.nii.gz', B1_map='B1.nii.gz',
                         noise=noise, verbose=vb).run()
        SteadyStateT1(sequence=seq, in_file=sim_file, verbose=vb).run()

        diff_M0 = Diff(in_file='SS_M0.nii.gz', baseline='M0.nii.gz',
                       noise = noise, verbose=vb).run()
        diff_T1 = Diff(in_file='SS_T1.nii.gz', baseline='T1.nii.gz',
                       noise=noise, verbose=vb).run()
        diff_B1 = Diff(in_file='SS_B1.nii.gz', baseline='B1.nii.gz',
                       noise=noise, verbose=vb).run()
        self.assertLessEqual(diff_M0.outputs.out_diff, 2)
        self.assertLessEqual(diff_T1.outputs.out_diff, 2)
        self.assertLessEqual(diff_B1.outputs.out_diff, 2)

if __name__ == '__main__':
    unittest.main()
