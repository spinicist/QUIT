from nipype.interfaces.base import CommandLineInputSpec, File, TraitedSpec, traits, isdefined
from nipype.interfaces.fsl.base import FSLCommand, FSLCommandInputSpec
from os import path


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
                     desc="folder of transformation matricies", xor=["tans_file"])
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
