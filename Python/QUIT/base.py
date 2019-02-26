# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""The quit module provides classes for interfacing with the `QUIT
<https://github.com/spinicist/QUIT>`_ command line tools.

Heavily inspired by the FSL base class
"""

import json
from os import path, getcwd
from nipype import logging
from nipype.utils.filemanip import fname_presuffix
from nipype.interfaces.base import traits, isdefined, CommandLine, CommandLineInputSpec, TraitedSpec, File, PackageInfo
from nipype.external.due import BibTeX


class InputBaseSpec(CommandLineInputSpec):
    """
    Base Input Specification for all QUIT Commands
    """

    # Inputs that are common to all program
    verbose = traits.Bool(desc='Print more information', argstr='--verbose')
    environ = {'QUIT_EXT': 'NIFTI_GZ'}


class InputSpec(InputBaseSpec):
    # Most QUIT programs take similar arguments

    # Most QUIT programs take a JSON parameter file
    param_file = File(desc='.json file', argstr='--json=%s',
                      xor=['param_dict'], mandatory=True, exists=True)

    param_dict = traits.Dict(desc='Param Dict', argstr='',
                             mandatory=True, xor=['param_file'])

    # Input nifti
    in_file = File(exists=True, argstr='%s', mandatory=True,
                   position=-1, desc='Input file')

    threads = traits.Int(
        desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(
        desc='Add a prefix to output filenames', argstr='--out=%s')
    mask_file = File(
        desc='Only process voxels within the mask', argstr='--mask=%s')
    residuals = traits.Bool(
        desc='Write out residuals for each data-point', argstr='--resids')


class BaseCommand(CommandLine):
    """
    Base support for QUIT commands.
    """

    references_ = [{
        'entry':
        BibTeX(
            '@article{Wood2017',
            'author = {Wood, Tobias Charles},'
            'doi = {https://doi.org/10.21105/joss.00656},'
            'journal = {Journal of Open Source Software},'
            'mendeley-groups = {Phd/Relaxometry},'
            'number = {26},'
            'pages = {656},'
            'title = {{QUIT: QUantitative Imaging Tools}},'
            'url = {https://github.com/spinicist/QUIT},'
            'volume = {3},'
            'year = {2017}'
            '}'),
        'tags': ['implementation'],
    }]

    def _make_json(self, value):
        fname = '_tmp_input.json'
        with open(fname, 'w') as outfile:
            json.dump(value, outfile)
        newarg = "--json=" + fname
        return newarg

    def _gen_fname(self,
                   basename,
                   cwd=None,
                   prefix=None,
                   suffix=None):
        """Generate a filename based on the given parameters.

        The filename will take the form: cwd/basename<suffix><ext>.
        If change_ext is True, it will use the extentions specified in
        <instance>intputs.output_type.

        Parameters
        ----------
        basename : str
            Filename to base the new filename on.
        cwd : str
            Path to prefix to the new filename. (default is os.getcwd())
        suffix : str
            Suffix to add to the `basename`.  (defaults is '' )
        change_ext : bool
            Flag to change the filename extension to the FSL output type.
            (default True)

        Returns
        -------
        fname : str
            New filename based on given parameters.

        """

        if basename == '':
            msg = 'Unable to generate filename for command %s. ' % self.cmd
            msg += 'basename is not set!'
            raise ValueError(msg)
        if cwd is None:
            cwd = getcwd()
        if prefix is None or not isdefined(prefix):
            prefix = ''
        if suffix is None:
            suffix = ''
        fname = fname_presuffix(
            basename, prefix=prefix, suffix=suffix, use_ext=True, newpath=cwd)
        return fname

    def _format_arg(self, name, spec, value):
        """
        Make parameter dictionary into a .json file for input to interface
        """
        if name == 'param_dict':
            return self._make_json(value)
        return super()._format_arg(name, spec, value)


class Command(BaseCommand):
    """
    Support for standard fitting tools.
    """

    def _add_prefixes(self, outputs):
        """
        Add prefix to everything in .outputs
        """
        for key, value in outputs.items():
            outputs[key] = self._gen_fname(value, prefix=self.inputs.prefix)
        return outputs

    def _list_outputs(self):
        return self._add_prefixes(self.output_spec().get())

####################### Simulator Base Classes #################################


class SimInputSpec(InputBaseSpec):
    # Most QUIT programs take similar arguments

    # Most QUIT programs take a JSON parameter file
    param_file = File(desc='.json file', argstr='--json=%s',
                      xor=['param_dict'], mandatory=True, exists=True)

    param_dict = traits.Dict(desc='Param Dict', argstr='',
                             mandatory=True, xor=['param_file'])

    # Input nifti
    in_file = File(exists=False, argstr='%s', mandatory=True,
                   position=-1, desc='Input file')

    threads = traits.Int(
        desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(
        desc='Add a prefix to output filenames', argstr='--out=%s')
    mask_file = File(
        desc='Only process voxels within the mask', argstr='--mask=%s')
    noise = traits.Float(desc='Noise level to add to simulation',
                         argstr='--simulate=%f', default_value=0.0, usedefault=True)


class SimOutputSpec(TraitedSpec):
    out_file = File(desc="Path to output image")


class SimCommand(Command):
    """
    Base support for QUIT simulator commands.
    """

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = path.abspath(self.inputs.in_file)
        return outputs


def check_QUIT():
    """
    Check if QUIT is installed
    """
    pass


def version_QUIT():
    """
    Check QUIT version
    """
    pass
