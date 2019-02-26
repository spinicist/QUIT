import unittest
from QUIT.base import BaseCommand
from QUIT.core import NewImage, Diff
from QUIT.relaxometry import DESPOT1, DESPOT1Sim, DESPOT2, DESPOT2Sim, HIFI, HIFISim, FM, FMSim

vb = True
BaseCommand.terminal_output = 'allatonce'


class DESPOT_SC(unittest.TestCase):
    def test_despot1(self):
        json = {'SPGR': {'TR': 10e-3, 'FA': [3, 18]},
                'T1File': 'T1.nii.gz', 'PDFile': 'PD.nii.gz'}
        spgr_file = 'sim_spgr.nii.gz'
        img_sz = [32, 32, 32]
        noise = 0.001

        NewImage(img_size=img_sz, grad_dim=0, grad_vals=(
            0.8, 1.0), out_file='PD.nii.gz', verbose=vb).run()
        NewImage(img_size=img_sz, grad_dim=1, grad_vals=(
            0.8, 1.3), out_file='T1.nii.gz', verbose=vb).run()

        DESPOT1Sim(param_dict=json, in_file=spgr_file,
                   noise=noise, verbose=vb).run()
        DESPOT1(param_dict=json, in_file=spgr_file, verbose=vb).run()

        diff_T1 = Diff(in_file='D1_T1.nii.gz', baseline='T1.nii.gz',
                       noise=noise, verbose=vb).run()
        diff_PD = Diff(in_file='D1_PD.nii.gz', baseline='PD.nii.gz',
                       noise=noise, verbose=vb).run()
        self.assertLessEqual(diff_T1.outputs.out_diff, 35)
        self.assertLessEqual(diff_PD.outputs.out_diff, 35)

    def test_hifi(self):
        json = {'SPGR': {'TR': 5e-3, 'FA': [3, 18]},
                'MPRAGE': {'FA': 5, 'TR': 5e-3, 'TI': 0.45, 'TD': 0, 'eta': 1, 'ETL': 64, 'k0': 0},
                'T1File': 'T1.nii.gz', 'PDFile': 'PD.nii.gz', 'B1File': 'B1.nii.gz'}
        spgr_file = 'sim_spgr.nii.gz'
        mprage_file = 'sim_mprage.nii.gz'
        img_sz = [32, 32, 32]
        noise = 0.001

        NewImage(out_file='PD.nii.gz', img_size=img_sz, grad_dim=0,
                 grad_vals=(0.8, 1.0), verbose=vb).run()
        NewImage(out_file='T1.nii.gz', img_size=img_sz, grad_dim=1,
                 grad_vals=(0.5, 1.5),  verbose=vb).run()
        NewImage(out_file='B1.nii.gz', img_size=img_sz, grad_dim=2,
                 grad_vals=(0.8, 1.2), verbose=vb).run()

        HIFISim(param_dict=json, spgr_file=spgr_file, mprage_file=mprage_file,
                noise=noise, verbose=vb).run()
        HIFI(param_dict=json, spgr_file=spgr_file,
             mprage_file=mprage_file, verbose=vb).run()

        diff_T1 = Diff(in_file='HIFI_T1.nii.gz', baseline='T1.nii.gz',
                       noise=noise, verbose=vb).run()
        diff_PD = Diff(in_file='HIFI_PD.nii.gz', baseline='PD.nii.gz',
                       noise=noise, verbose=vb).run()
        diff_B1 = Diff(in_file='HIFI_B1.nii.gz', baseline='B1.nii.gz',
                       noise=noise, verbose=vb).run()
        self.assertLessEqual(diff_T1.outputs.out_diff, 60)
        self.assertLessEqual(diff_PD.outputs.out_diff, 35)
        self.assertLessEqual(diff_B1.outputs.out_diff, 60)

    def test_despot2(self, gs=False, tol=20):
        json = {'SSFP': {'TR': 10e-3, 'FA': [15, 30, 45, 60], 'PhaseInc': [
            180, 180, 180, 180]}, 'T2File': 'T2.nii.gz', 'PDFile': 'PD.nii.gz'}
        ssfp_file = 'sim_ssfp.nii.gz'
        img_sz = [32, 32, 32]
        noise = 0.001

        NewImage(img_size=img_sz, grad_dim=0, grad_vals=(0.8, 1.0),
                 out_file='PD.nii.gz', verbose=vb).run()
        NewImage(img_size=img_sz, grad_dim=1, grad_vals=(0.8, 1.2),
                 out_file='T1.nii.gz', verbose=vb).run()
        NewImage(img_size=img_sz, grad_dim=2, grad_vals=(0.04, 0.1),
                 out_file='T2.nii.gz', verbose=vb).run()

        DESPOT2Sim(param_dict=json, in_file=ssfp_file,
                   t1_file='T1.nii.gz', ellipse=gs, noise=noise, verbose=vb).run()
        DESPOT2(param_dict=json, in_file=ssfp_file,
                t1_file='T1.nii.gz', ellipse=gs, verbose=vb).run()

        diff_T2 = Diff(in_file='D2_T2.nii.gz', baseline='T2.nii.gz',
                       noise=noise, verbose=vb).run()
        diff_PD = Diff(in_file='D2_PD.nii.gz', baseline='PD.nii.gz',
                       noise=noise, verbose=vb).run()
        self.assertLessEqual(diff_T2.outputs.out_diff, tol)
        self.assertLessEqual(diff_PD.outputs.out_diff, tol)

    def test_despot2gs(self):
        self.test_despot2(True, 30)

    def test_fm(self):
        json = {'SSFP': {'TR': 5e-3, 'FA': [15, 15, 60, 60], 'PhaseInc': [180, 0, 180, 0]},
                'T2File': 'T2.nii.gz', 'PDFile': 'PD.nii.gz', 'f0File': 'f0.nii.gz'}
        ssfp_file = 'sim_ssfp.nii.gz'
        img_sz = [16, 16, 16]
        noise = 0.001

        NewImage(img_size=img_sz, fill=1.0,
                 out_file='PD.nii.gz', verbose=vb).run()
        NewImage(img_size=img_sz, grad_dim=0, grad_vals=(0.8, 1.2),
                 out_file='T1.nii.gz', verbose=vb).run()
        NewImage(img_size=img_sz, grad_dim=1, grad_vals=(0.04, 0.1),
                 out_file='T2.nii.gz', verbose=vb).run()
        NewImage(img_size=img_sz, grad_dim=2, grad_vals=(-100, 100),
                 out_file='f0.nii.gz', verbose=vb).run()

        FMSim(param_dict=json, in_file=ssfp_file, asym=False,
              t1_file='T1.nii.gz', noise=noise, verbose=vb).run()
        FM(param_dict=json, in_file=ssfp_file, asym=False,
           t1_file='T1.nii.gz', verbose=vb).run()

        diff_T2 = Diff(in_file='FM_T2.nii.gz', baseline='T2.nii.gz',
                       noise=noise, verbose=vb).run()
        diff_PD = Diff(in_file='FM_PD.nii.gz', baseline='PD.nii.gz',
                       noise=noise, verbose=vb).run()
        self.assertLessEqual(diff_T2.outputs.out_diff, 20)
        self.assertLessEqual(diff_PD.outputs.out_diff, 10)


if __name__ == '__main__':
    unittest.main()
