import unittest
from nipype.interfaces.base import CommandLine
from QUIT.core import NewImage, Diff
from QUIT.relaxometry import Multiecho, MultiechoSim, MUPA, MUPASim

vb = True
CommandLine.terminal_output = 'allatonce'


class Relax(unittest.TestCase):
    def test_multiecho(self):
        me = {'MultiEcho': {'TR': 10, 'TE1': 0.01, 'ESP': 0.01,
                            'ETL': 5}}
        me_file = 'sim_me.nii.gz'
        img_sz = [32, 32, 32]
        noise = 0.001

        NewImage(img_size=img_sz, grad_dim=0, grad_vals=(0.8, 1.0),
                 out_file='PD.nii.gz', verbose=vb).run()
        NewImage(img_size=img_sz, grad_dim=2, grad_vals=(0.04, 0.1),
                 out_file='T2.nii.gz', verbose=vb).run()

        MultiechoSim(sequence=me, in_file=me_file,
                     PD='PD.nii.gz', T2='T2.nii.gz',
                     noise=noise, verbose=vb).run()
        Multiecho(sequence=me, in_file=me_file, verbose=vb).run()

        diff_T2 = Diff(in_file='ME_T2.nii.gz', baseline='T2.nii.gz',
                       noise=noise, verbose=vb).run()
        diff_PD = Diff(in_file='ME_PD.nii.gz', baseline='PD.nii.gz',
                       noise=noise, verbose=vb).run()
        self.assertLessEqual(diff_T2.outputs.out_diff, 3)
        self.assertLessEqual(diff_PD.outputs.out_diff, 2)

    def test_MUPA(self):
        seq = {'MUPA': {'TR': 0.002786, 'Tramp': 0.01, 'SPS': 64, 'FA': 2,
                        "prep_type": ["inversion", "delay", "delay", "delay", "t2-prep", "t2-prep", "t2-prep"],
                        "prep_time": [0.03, 0, 0, 0, 0.04, 0.08, 0.16]}}
        sim_file = 'sim_mupa.nii.gz'
        img_sz = [16, 16, 16]
        noise = 0.01

        NewImage(img_size=img_sz, grad_dim=0, grad_vals=(80, 120),
                 out_file='PD.nii.gz', verbose=vb).run()
        NewImage(img_size=img_sz, grad_dim=0, grad_vals=(0.8, 3.0),
                 out_file='T1.nii.gz', verbose=vb).run()
        NewImage(img_size=img_sz, grad_dim=0, grad_vals=(0.04, 1.0),
                 out_file='T2.nii.gz', verbose=vb).run()

        MUPASim(sequence=seq, in_file=sim_file,
                PD='PD.nii.gz', T1='T1.nii.gz', T2='T2.nii.gz',
                noise=noise, verbose=vb).run()
        MUPA(sequence=seq, in_file=sim_file, verbose=vb, threads=-1).run()

        diff_T2 = Diff(in_file='MUPA_T2.nii.gz', baseline='T2.nii.gz',
                       noise=noise, verbose=vb).run()
        diff_PD = Diff(in_file='MUPA_PD.nii.gz', baseline='PD.nii.gz',
                       noise=noise, verbose=vb).run()
        diff_T1 = Diff(in_file='MUPA_T1.nii.gz', baseline='T1.nii.gz',
                       noise=noise, verbose=vb).run()

        self.assertLessEqual(diff_PD.outputs.out_diff, 10)
        self.assertLessEqual(diff_T1.outputs.out_diff, 10)
        self.assertLessEqual(diff_T2.outputs.out_diff, 12)


if __name__ == '__main__':
    unittest.main()
