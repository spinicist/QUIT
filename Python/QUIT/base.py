# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""The quit module provides classes for interfacing with the `QUIT
<https://github.com/spinicist/QUIT>`_ command line tools.

Heavily inspired by the FSL base class
"""

import json
from copy import deepcopy
from os import path, getcwd
from nipype import logging
from nipype.utils.filemanip import fname_presuffix
from nipype.interfaces.base import traits, isdefined, CommandLine, CommandLineInputSpec, TraitedSpec, DynamicTraitedSpec, File, PackageInfo
from nipype.external.due import BibTeX

################################# Input Specs ##################################


class InputBaseSpec(CommandLineInputSpec):
    """
    Barebones input spec suitable for utilities like NewImage
    """
    verbose = traits.Bool(desc='Print more information', argstr='--verbose')
    environ = {'QUIT_EXT': 'NIFTI_GZ'}


class InputSpec(InputBaseSpec):
    """
    Input Specification for most QUIT tools
    """
    # Most QUIT programs take similar arguments
    json = traits.File(exists=True, desc='JSON Input file', argstr='--json=%s')

    subregion = traits.String(
        desc='Only process a subregion of the image. Argument should be a string "start_x,start_y,start_z,size_x,size_y,size_z"', argstr='--subregion=%s')
    threads = traits.Int(
        desc='Use N threads (default=4, 0=hardware limit)', argstr='--threads=%d')
    prefix = traits.String(
        desc='Add a prefix to output filenames', argstr='--out=%s')
    mask_file = File(
        desc='Only process voxels within the mask', argstr='--mask=%s')


class FitInputSpec(InputSpec):
    """
    Input Specification for tools in fitting mode
    """
    # Input nifti
    in_file = File(exists=True, argstr='%s', mandatory=True,
                   position=-1, desc='Input file')
    covar = traits.Bool(
        decsc='Write out parameter covar images', argstr='--covar')
    residuals = traits.Bool(
        desc='Write out residuals for each data-point', argstr='--resids')


class SimInputBaseSpec(DynamicTraitedSpec):
    """
    Input specification for tools in simulation mode
    """
    ignore_exception = traits.Bool(
        False,
        usedefault=True,
        nohash=True,
        deprecated='1.0.0',
        desc='Print an error message instead of throwing an exception '
        'in case the interface fails to run')
    args = traits.Str(argstr='%s', desc='Additional parameters to the command')
    environ = traits.DictStrStr(
        desc='Environment variables', usedefault=True, nohash=True)
    # This input does not have a "usedefault=True" so the set_default_terminal_output()
    # method would work
    terminal_output = traits.Enum(
        'stream',
        'allatonce',
        'file',
        'none',
        deprecated='1.0.0',
        desc=('Control terminal output: `stream` - '
              'displays to terminal immediately (default), '
              '`allatonce` - waits till command is '
              'finished to display output, `file` - '
              'writes output to file, `none` - output'
              ' is ignored'),
        nohash=True)

    json = traits.File(exists=True, desc='JSON Input file', argstr='--json=%s')
    noise = traits.Float(desc='Noise level to add to simulation',
                         argstr='--simulate=%f', default_value=0.0, usedefault=True)


class SimInputSpec(SimInputBaseSpec):
    """
    Input specification for tools in simulation mode
    """
    in_file = File(argstr='%s', mandatory=True,
                   position=-1, desc='Output simulated file')
################################# Output Specs #################################


class SimOutputSpec(TraitedSpec):
    out_file = File(desc="Path to output image")

################################### Commands ###################################


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

    def _parse_inputs(self, skip=None):
        # Make sequence dictionary into a .json file for input to interface
        if hasattr(self, '_json') and not isdefined(self.inputs.json):
            fname = '_input.json'
            with open(fname, 'w') as outfile:
                json.dump(self._json, outfile, indent=2)
            self.inputs.json = path.abspath(fname)
        return super()._parse_inputs(skip)


class FitCommand(BaseCommand):
    """
    Support for standard fitting tools.
    """

    def _add_prefixes(self, outputs):
        """
        Add prefix to everything in .outputs
        """
        for key, value in outputs.items():
            if isinstance(value, (list,)):
                outputs[key] = [self._gen_fname(
                    x, prefix=self.inputs.prefix) for x in value]
            else:
                outputs[key] = self._gen_fname(
                    value, prefix=self.inputs.prefix)
        return outputs

    def _list_outputs(self, outputs=None):
        if outputs is None:
            outputs = self.output_spec().get()
        return self._add_prefixes(outputs)

    def __init__(self, sequence={}, **kwargs):
        self._json = deepcopy(sequence)
        super().__init__(**kwargs)


class SimCommand(BaseCommand):
    """
    Base support for QUIT simulator commands.
    """

    def __init__(self, sequence={}, **kwargs):
        self._json = deepcopy(sequence)
        for pf in self._param_files:
            setattr(self.input_spec, pf, File(
                exists=True, desc='Path to %s map' % pf))
        super().__init__(**kwargs)

    def _parse_inputs(self, skip=[]):
        for pf in self._param_files:
            if isdefined(getattr(self.inputs, pf)):
                self._json[pf] = getattr(self.inputs, pf)
                skip.append(pf)
        return super()._parse_inputs(skip)

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
