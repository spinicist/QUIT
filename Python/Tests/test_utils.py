from pathlib import Path
from os import chdir
import unittest
from math import sqrt
from nipype.interfaces.base import CommandLine
from qipype.interfaces.core import NewImage, Diff
from qipype.interfaces.utils import PolyImage, PolyFit, Filter, RFProfile

vb = True
CommandLine.terminal_output = 'allatonce'


class Utils(unittest.TestCase):
    def setUp(self):
        Path('testdata').mkdir(exist_ok=True)
        chdir('testdata')

    def tearDown(self):
        chdir('../')

    def test_poly(self):
        sz = 16
        scale = sqrt(3 * (sz/2)**2)
        poly = {'center': [0, 0, 0],
                'scale':  scale,
                'coeffs': [1, 1, 2, 4, 1, 0, 0, 1, 0, 1]
                }
        poly_sim = 'poly_sim.nii.gz'
        mask_file = 'poly_mask.nii.gz'

        NewImage(img_size=[sz, sz, sz], grad_dim=0, grad_vals=(0, 1), grad_steps=1,
                 out_file=mask_file, verbose=vb).run()
        PolyImage(ref_file=mask_file, out_file=poly_sim,
                  order=2, poly=poly, verbose=vb).run()
        fit = PolyFit(in_file=poly_sim, order=2, robust=True).run()

        terms_diff = sum([abs(x - y) for x,
                          y in zip(poly['coeffs'], fit.outputs.poly['coeffs'])])
        self.assertLessEqual(terms_diff, 1.e-6)

        PolyImage(ref_file=mask_file, out_file='poly_sim2.nii.gz',
                  order=2, poly=fit.outputs.poly, verbose=vb).run()
        img_diff = Diff(baseline=poly_sim,
                        in_file='poly_sim2.nii.gz', noise=1).run()
        self.assertLessEqual(img_diff.outputs.out_diff, 1.e-3)

    def test_kfilter(self):
        NewImage(out_file='steps.nii.gz', img_size=[64, 64, 64],
                 grad_dim=0, grad_vals=(0, 8), grad_steps=4, verbose=vb).run()
        Filter(in_file='steps.nii.gz', filter_spec='Gauss,2.0', verbose=vb).run()

    def test_rfprofile(self):
        NewImage(out_file='rf_b1plus.nii.gz', img_size=[32, 32, 32],
                 fill=1.0, verbose=vb).run()
        RFProfile(in_file='rf_b1plus.nii.gz', out_file='rf_slab.nii.gz',
                  rf={'rf_pos': [0, 1], 'rf_vals': [0, 1]}, verbose=vb).run()
        NewImage(out_file='rf_ref.nii.gz', img_size=[32, 32, 32], grad_dim=2,
                 grad_vals=(-16, 15), verbose=vb).run()
        rf_diff = Diff(baseline='rf_ref.nii.gz', in_file='rf_slab.nii.gz',
                       noise=1, abs_diff=True, verbose=vb).run()
        self.assertLessEqual(rf_diff.outputs.out_diff, 1.e-3)


if __name__ == '__main__':
    unittest.main()
