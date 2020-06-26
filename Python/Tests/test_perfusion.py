import unittest
from math import sqrt
from nipype.interfaces.base import CommandLine
from qipype.interfaces.core import NewImage, Diff
from qipype.interfaces.perfusion import ASL, ASE, ASESim, ZShim

vb = True
CommandLine.terminal_output = 'allatonce'


class Perfusion(unittest.TestCase):
    def test_asl(self):
        seq = {'CASL': {'TR': 4.0, 'label_time': 3.0,
                        'post_label_delay': [0.3]}}
        asl_file = 'sim_asl.nii.gz'

        NewImage(img_size=[32, 32, 32, 2], grad_dim=3, grad_vals=(1, 1.06), grad_steps=1,
                 out_file=asl_file, verbose=vb).run()
        NewImage(img_size=[32, 32, 32], fill=147.355,
                 out_file='ref_cbf.nii.gz', verbose=vb).run()

        ASL(sequence=seq, in_file=asl_file, verbose=vb).run()

        diff_CBF = Diff(in_file='CASL_CBF.nii.gz',
                        baseline='ref_cbf.nii.gz', verbose=vb).run()
        self.assertLessEqual(diff_CBF.outputs.out_diff, 0.1)

    def test_oef(self):
        # Use MultiEchoFlex as a proxy for ASE
        seq = {'MultiEcho': {'TR': 2.5,
                             'TE': [-0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07]}}
        ase_file = 'sim_ase.nii.gz'
        img_sz = [32, 32, 32]
        noise = 0.001

        NewImage(img_size=img_sz, grad_dim=0, fill=100.,
                 out_file='S0.nii.gz', verbose=vb).run()
        NewImage(img_size=img_sz, grad_dim=0, grad_vals=(-0.01, 0.01),
                 out_file='dT.nii.gz', verbose=vb).run()
        NewImage(img_size=img_sz, grad_dim=1, grad_vals=(1.0, 3.0),
                 out_file='R2p.nii.gz', verbose=vb).run()
        NewImage(img_size=img_sz, grad_dim=2, grad_vals=(0.005, 0.025),
                 out_file='DBV.nii.gz', verbose=vb).run()

        ASESim(sequence=seq, out_file=ase_file, noise=noise, verbose=vb,
               S0_map='S0.nii.gz',
               dT_map='dT.nii.gz',
               R2p_map='R2p.nii.gz',
               DBV_map='DBV.nii.gz').run()
        ASE(sequence=seq, in_file=ase_file, verbose=vb).run()

        diff_R2p = Diff(in_file='ASE_R2p.nii.gz', baseline='R2p.nii.gz',
                        noise=noise, verbose=vb).run()
        diff_DBV = Diff(in_file='ASE_DBV.nii.gz', baseline='DBV.nii.gz',
                        noise=noise, verbose=vb).run()
        self.assertLessEqual(diff_R2p.outputs.out_diff, 12)
        self.assertLessEqual(diff_DBV.outputs.out_diff, 50)

    def test_oef_fixed_dbv(self):
        # Use MultiEchoFlex as a proxy for ASE
        seq = {'MultiEcho': {'TR': 2.5,
                             'TE': [-0.05, -0.04, -0.03, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.04, 0.05]}}
        ase_file = 'sim_ase.nii.gz'
        img_sz = [32, 32, 32]
        noise = 0.001
        DBV = 0.01

        NewImage(img_size=img_sz, grad_dim=0, fill=100.,
                 out_file='S0.nii.gz', verbose=vb).run()
        NewImage(img_size=img_sz, grad_dim=0, grad_vals=(-0.01, 0.01),
                 out_file='dT.nii.gz', verbose=vb).run()
        NewImage(img_size=img_sz, grad_dim=1, grad_vals=(1.0, 3.0),
                 out_file='R2p.nii.gz', verbose=vb).run()

        ASESim(sequence=seq, out_file=ase_file,
               fix_DBV=DBV, noise=noise, verbose=vb,
               S0_map='S0.nii.gz',
               dT_map='dT.nii.gz',
               R2p_map='R2p.nii.gz').run()
        ASE(sequence=seq, in_file=ase_file, fix_DBV=DBV, verbose=vb).run()

        diff_R2p = Diff(in_file='ASE_R2p.nii.gz', baseline='R2p.nii.gz',
                        noise=noise, verbose=vb).run()
        self.assertLessEqual(diff_R2p.outputs.out_diff, 1.0)

    def test_zshim(self):
        nshims = 8
        sz = 32
        ref_val = sqrt(sum([x**2 for x in range(1, nshims + 1)]))
        NewImage(out_file='zshim.nii.gz', img_size=[sz, sz, sz, nshims], grad_dim=3,
                 grad_vals=(1, nshims), grad_steps=7,  verbose=vb).run()
        NewImage(out_file='zshim_ref.nii.gz', img_size=[sz, sz, sz],
                 fill=ref_val,  verbose=vb).run()
        ZShim(in_file='zshim.nii.gz', zshims=nshims, verbose=vb).run()
        zdiff = Diff(in_file='zshim_zshim.nii.gz',
                     baseline='zshim_ref.nii.gz', noise=1, verbose=vb).run()
        self.assertLessEqual(zdiff.outputs.out_diff, 1.e-3)


if __name__ == '__main__':
    unittest.main()
