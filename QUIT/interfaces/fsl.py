from nipype.interfaces.base import CommandLineInputSpec, Directory, File, TraitedSpec, OutputMultiPath, traits, isdefined
from nipype.interfaces.fsl.base import FSLCommand, FSLCommandInputSpec, Info
from nipype import LooseVersion
from nibabel import load
import os
from os import path


class MCFLIRTInputSpec(FSLCommandInputSpec):
    in_file = File(
        exists=True,
        position=0,
        argstr="-in %s",
        mandatory=True,
        desc="timeseries to motion-correct")
    out_file = File(
        argstr='-out %s', genfile=True, desc="file to write", hash_files=False)
    cost = traits.Enum(
        'mutualinfo',
        'woods',
        'corratio',
        'normcorr',
        'normmi',
        'leastsquares',
        argstr='-cost %s',
        desc="cost function to optimize")
    bins = traits.Int(argstr='-bins %d', desc="number of histogram bins")
    dof = traits.Int(
        argstr='-dof %d', desc="degrees of freedom for the transformation")
    ref_vol = traits.Int(argstr='-refvol %d', desc="volume to align frames to")
    scaling = traits.Float(
        argstr='-scaling %.2f', desc="scaling factor to use")
    smooth = traits.Float(
        argstr='-smooth %.2f', desc="smoothing factor for the cost function")
    rotation = traits.Int(
        argstr='-rotation %d', desc="scaling factor for rotation tolerances")
    stages = traits.Int(
        argstr='-stages %d',
        desc="stages (if 4, perform final search with sinc interpolation")
    init = File(
        exists=True, argstr='-init %s', desc="inital transformation matrix")
    interpolation = traits.Enum(
        "spline",
        "nn",
        "sinc",
        argstr="-%s_final",
        desc="interpolation method for transformation")
    use_gradient = traits.Bool(
        argstr='-gdt', desc="run search on gradient images")
    use_contour = traits.Bool(
        argstr='-edge', desc="run search on contour images")
    mean_vol = traits.Bool(argstr='-meanvol', desc="register to mean volume")
    stats_imgs = traits.Bool(
        argstr='-stats', desc="produce variance and std. dev. images")
    save_mats = traits.Bool(
        argstr='-mats', desc="save transformation matrices")
    save_plots = traits.Bool(
        argstr='-plots', desc="save transformation parameters")
    save_rms = traits.Bool(
        argstr='-rmsabs -rmsrel', desc="save rms displacement parameters")
    ref_file = File(
        exists=True,
        argstr='-reffile %s',
        desc="target image for motion correction")


class MCFLIRTOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="motion-corrected timeseries")
    variance_img = File(exists=True, desc="variance image")
    std_img = File(exists=True, desc="standard deviation image")
    mean_img = File(
        exists=True, desc="mean timeseries image (if mean_vol=True)")
    par_file = File(exists=True, desc="text-file with motion parameters")
    mat_dir = Directory(
        exists=True, desc="directory containing transformation matrices")
    mat_file = OutputMultiPath(
        File(exists=True), desc="transformation matrices")
    rms_files = OutputMultiPath(
        File(exists=True),
        desc="absolute and relative displacement parameters")


class MCFLIRT(FSLCommand):
    """FSL MCFLIRT wrapper for within-modality motion correction

    For complete details, see the `MCFLIRT Documentation.
    <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/MCFLIRT>`_

    Examples
    --------
    >>> from nipype.interfaces import fsl
    >>> mcflt = fsl.MCFLIRT()
    >>> mcflt.inputs.in_file = 'functional.nii'
    >>> mcflt.inputs.cost = 'mutualinfo'
    >>> mcflt.inputs.out_file = 'moco.nii'
    >>> mcflt.cmdline
    'mcflirt -in functional.nii -cost mutualinfo -out moco.nii'
    >>> res = mcflt.run()  # doctest: +SKIP

    """
    _cmd = 'mcflirt'
    input_spec = MCFLIRTInputSpec
    output_spec = MCFLIRTOutputSpec

    def _format_arg(self, name, spec, value):
        if name == "interpolation":
            if value == "trilinear":
                return ""
            else:
                return spec.argstr % value
        return super(MCFLIRT, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self._outputs().get()

        outputs['out_file'] = self._gen_outfilename()
        output_dir = os.path.dirname(outputs['out_file'])

        if isdefined(self.inputs.stats_imgs) and self.inputs.stats_imgs:
            if LooseVersion(Info.version()) < LooseVersion('6.0.0'):
                # FSL <6.0 outputs have .nii.gz_variance.nii.gz as extension
                outputs['variance_img'] = self._gen_fname(
                    outputs['out_file'] + '_variance.ext', cwd=output_dir)
                outputs['std_img'] = self._gen_fname(
                    outputs['out_file'] + '_sigma.ext', cwd=output_dir)
            else:
                outputs['variance_img'] = self._gen_fname(
                    outputs['out_file'], suffix='_variance', cwd=output_dir)
                outputs['std_img'] = self._gen_fname(
                    outputs['out_file'], suffix='_sigma', cwd=output_dir)

        # The mean image created if -stats option is specified ('meanvol')
        # is missing the top and bottom slices. Therefore we only expose the
        # mean image created by -meanvol option ('mean_reg') which isn't
        # corrupted.
        # Note that the same problem holds for the std and variance image.

        if isdefined(self.inputs.mean_vol) and self.inputs.mean_vol:
            if LooseVersion(Info.version()) < LooseVersion('6.0.0'):
                # FSL <6.0 outputs have .nii.gz_mean_img.nii.gz as extension
                outputs['mean_img'] = self._gen_fname(
                    outputs['out_file'] + '_mean_reg.ext', cwd=output_dir)
            else:
                outputs['mean_img'] = self._gen_fname(
                    outputs['out_file'], suffix='_mean_reg', cwd=output_dir)

        if isdefined(self.inputs.save_mats) and self.inputs.save_mats:
            _, filename = os.path.split(outputs['out_file'])
            mat_dir = os.path.join(output_dir, filename + '.mat')
            outputs['mat_dir'] = mat_dir
            _, _, _, timepoints = load(self.inputs.in_file).shape
            outputs['mat_file'] = []
            for t in range(timepoints):
                outputs['mat_file'].append(
                    os.path.join(mat_dir, 'MAT_%04d' % t))
        if isdefined(self.inputs.save_plots) and self.inputs.save_plots:
            # Note - if e.g. out_file has .nii.gz, you get .nii.gz.par,
            # which is what mcflirt does!
            outputs['par_file'] = outputs['out_file'] + '.par'
        if isdefined(self.inputs.save_rms) and self.inputs.save_rms:
            outfile = outputs['out_file']
            outputs['rms_files'] = [outfile + '_abs.rms', outfile + '_rel.rms']
        return outputs

    def _gen_filename(self, name):
        if name == 'out_file':
            return self._gen_outfilename()
        return None

    def _gen_outfilename(self):
        out_file = self.inputs.out_file
        if isdefined(out_file):
            out_file = os.path.realpath(out_file)
        if not isdefined(out_file) and isdefined(self.inputs.in_file):
            out_file = self._gen_fname(self.inputs.in_file, suffix='_mcf')
        return os.path.abspath(out_file)


class ApplyXfm4DInputSpec(FSLCommandInputSpec):
    in_file = File(exists=True, position=0, argstr='%s',
                   mandatory=True, desc="timeseries to motion-correct")
    ref_vol = File(exists=True, position=1, argstr='%s',
                   mandatory=True, desc="volume with final FOV and resolution")
    out_file = File(exists=True, position=2, argstr='%s',
                    genfile=True, desc="file to write", hash_files=False)
    trans_file = File(argstr='%s', position=3, desc="single tranformation matrix", xor=[
        "trans_dir"], requires=["single_matrix"])
    trans_dir = File(argstr='%s', position=3,
                     desc="folder of transformation matrices", xor=["trans_file"])
    single_matrix = traits.Bool(
        argstr='-singlematrix', desc="true if applying one volume to all timepoints")
    four_digit = traits.Bool(
        argstr='-fourdigit', desc="true mat names have four digits not five")
    user_prefix = traits.Str(
        argstr='-userprefix %s', desc="supplied prefix if mats don't start with 'MAT_'")


class ApplyXfm4DOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="transform applied timeseries")


class ApplyXfm4D(FSLCommand):
    """
    Wraps the applyxfm4D command line tool for applying one 3D transform to every volume in a 4D image OR
    a directory of 3D tansforms to a 4D image of the same length.

    Examples
    ---------
    >>> import nipype.interfaces.fsl as fsl
    >>> from nipype.testing import example_data
    >>> applyxfm4d = fsl.ApplyXfm4D()
    >>> applyxfm4d.inputs.in_file = example_data('functional.nii')
    >>> applyxfm4d.inputs.in_matrix_file = example_data('functional_mcf.mat')
    >>> applyxfm4d.inputs.out_file = 'newfile.nii'
    >>> applyxfm4d.inputs.reference = example_data('functional_mcf.nii')
    >>> result = applyxfm.run() # doctest: +SKIP

    """

    _cmd = 'applyxfm4D'
    input_spec = ApplyXfm4DInputSpec
    output_spec = ApplyXfm4DOutputSpec

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_file'] = self._gen_outfilename()
        return outputs

    def _gen_filename(self, name):
        if name == 'out_file':
            return self._gen_outfilename()
        return None

    def _gen_outfilename(self):
        out_file = self.inputs.out_file
        if isdefined(out_file):
            out_file = path.realpath(out_file)
        if not isdefined(out_file) and isdefined(self.inputs.in_file):
            out_file = self._gen_fname(
                self.inputs.in_file, suffix='_warp4D')
        return path.abspath(out_file)
