import unittest
import numpy as np
from nipype.interfaces.base import CommandLine
from QUIT.core import NewImage, Diff
from QUIT.ssfp import Ellipse, EllipseSim, PLANET, PLANETSim, eMT, eMTSim

vb = True
CommandLine.terminal_output = 'allatonce'


class SSFP(unittest.TestCase):
    def test_planet(self):
        ellipse_seq = {"SSFP": {
            "FA": [15, 15, 15, 15, 15, 15],
            "PhaseInc": [180, 240, 300, 0, 60, 120],
            "TR": 0.01
        }
        }
        planet_seq = {
            "SSFP": {
                "FA": [15],
                "PhaseInc": [0],
                "TR": 0.01
            }
        }
        ellipse_file = 'planet_ellipse.nii.gz'
        planet_G = 'planet_G.nii.gz'
        planet_a = 'planet_a.nii.gz'
        planet_b = 'planet_b.nii.gz'
        img_sz = [32, 32, 32]
        noise = 0.001

        NewImage(out_file='PD.nii.gz', verbose=vb, img_size=img_sz,
                 fill=1).run()
        NewImage(out_file='T1.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=0, grad_vals=(0.8, 1.3)).run()
        NewImage(out_file='T2.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=1, grad_vals=(0.05, 0.1)).run()
        NewImage(out_file='zero.nii.gz', verbose=vb, img_size=img_sz,
                 fill=0).run()

        PLANETSim(sequence=planet_seq, G_file=planet_G, a_file=planet_a, b_file=planet_b,
                  noise=0, verbose=vb,
                  PD='PD.nii.gz',
                  T1='T1.nii.gz',
                  T2='T2.nii.gz').run()
        EllipseSim(sequence=ellipse_seq, in_file=ellipse_file,
                   noise=noise, verbose=vb,
                   G=planet_G, a=planet_a, b=planet_b, theta_0='zero.nii.gz', phi_rf='zero.nii.gz').run()
        Ellipse(sequence=ellipse_seq, in_file=ellipse_file, verbose=vb).run()
        PLANET(sequence=planet_seq, G_map=planet_G,
               a_map=planet_a, b_map=planet_b, verbose=vb).run()

        diff_G = Diff(in_file='ES_G.nii.gz', baseline=planet_G,
                      noise=noise, verbose=vb).run()
        diff_a = Diff(in_file='ES_a.nii.gz', baseline=planet_a,
                      noise=noise, verbose=vb).run()
        diff_b = Diff(in_file='ES_b.nii.gz', baseline=planet_b,
                      noise=noise, verbose=vb).run()
        diff_PD = Diff(in_file='PLANET_PD.nii.gz', baseline='PD.nii.gz',
                       noise=noise, verbose=vb).run()
        diff_T1 = Diff(in_file='PLANET_T1.nii.gz', baseline='T1.nii.gz',
                       noise=noise, verbose=vb).run()
        diff_T2 = Diff(in_file='PLANET_T2.nii.gz', baseline='T2.nii.gz',
                       noise=noise, verbose=vb).run()
        self.assertLessEqual(diff_G.outputs.out_diff, 3.5)
        self.assertLessEqual(diff_a.outputs.out_diff, 3)
        self.assertLessEqual(diff_b.outputs.out_diff, 7)
        self.assertLessEqual(diff_PD.outputs.out_diff, 4)
        self.assertLessEqual(diff_T1.outputs.out_diff, 1)
        self.assertLessEqual(diff_T2.outputs.out_diff, 1)

    def test_emt(self):
        ellipse_sim = {"SSFP": {
            "FA": [1, 1, 1, 1, 1, 1],
            "PhaseInc": [180, 240, 300, 0, 60, 120],
            "TR": 0.01
        }}
        ellipse_fit = {"SSFP": {
            "FA": [1, 1, 1, 1, 1, 1],
            "PhaseInc": [180, 240, 300, 0, 60, 120],
            "TR": 0.01
        }}
        emt_seq = {
            "SSFPMT": {
                "FA": [30],
                "PhaseInc": [180, 240, 300, 0, 60, 120],
                "TR": [0.01],
                "Trf": [0.00025],
                "pulse": {"p1": 0.250629, "p2": 0.201183, "bandwidth": 2}
            }
        }
        ellipse_file = 'emt_ellipse.nii.gz'
        emt_G = 'emt_G.nii.gz'
        emt_a = 'emt_a.nii.gz'
        emt_b = 'emt_b.nii.gz'
        img_sz = [32, 32, 32]
        noise = 0.001

        NewImage(out_file='PD.nii.gz', verbose=vb, img_size=img_sz,
                 fill=1).run()
        NewImage(out_file='T1_f.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=0, grad_vals=(0.8, 1.3)).run()
        NewImage(out_file='T2_f.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=1, grad_vals=(0.05, 0.1)).run()
        NewImage(out_file='f_b.nii.gz', verbose=vb, img_size=img_sz,
                 grad_dim=2, grad_vals=(0.01, 0.15)).run()
        NewImage(out_file='k_bf.nii.gz', verbose=vb, img_size=img_sz,
                 fill=2).run()
        NewImage(out_file='zero.nii.gz', verbose=vb, img_size=img_sz,
                 fill=0).run()

        eMTSim(sequence=emt_seq, G_file=emt_G, a_file=emt_a, b_file=emt_b,
               noise=0, verbose=vb,
               PD='PD.nii.gz',
               T1_f='T1_f.nii.gz',
               T2_f='T2_f.nii.gz',
               f_b='f_b.nii.gz',
               k_bf='k_bf.nii.gz').run()
        EllipseSim(sequence=ellipse_sim, in_file=ellipse_file,
                   noise=noise, verbose=vb,
                   G=emt_G, a=emt_a, b=emt_b, theta_0='zero.nii.gz', phi_rf='zero.nii.gz').run()
        Ellipse(sequence=ellipse_fit, in_file=ellipse_file, verbose=vb).run()
        eMT(sequence=emt_seq, G_map=emt_G,
            a_map=emt_a, b_map=emt_b, verbose=vb).run()

        # Currently the simulation framework does not support blocked algorithms
        # Hence simulating a proper eMT dataset would involve individually simulating
        # the different TR/Trf combinations and then merging them. I'm not doing that
        # now, so a proper test will have to wait

        diff_G = Diff(in_file='ES_G.nii.gz', baseline=emt_G,
                      noise=noise, verbose=vb).run()
        diff_a = Diff(in_file='ES_a.nii.gz', baseline=emt_a,
                      noise=noise, verbose=vb).run()
        diff_b = Diff(in_file='ES_b.nii.gz', baseline=emt_b,
                      noise=noise, verbose=vb).run()
        self.assertLessEqual(diff_G.outputs.out_diff, 3.5)
        self.assertLessEqual(diff_a.outputs.out_diff, 3.5)
        self.assertLessEqual(diff_b.outputs.out_diff, 10)


if __name__ == '__main__':
    unittest.main()
