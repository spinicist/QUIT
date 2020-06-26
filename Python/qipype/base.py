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


def SimInputSpec(name, varying, fixed=[], out_files=None, extras=None):
    """
    Input specification for tools in simulation mode
    """
    if out_files is None:
        out_files = ['out']

    attrs = {'noise': traits.Float(desc='Noise level to add to simulation',
                                   argstr='--simulate=%f', default_value=0.0, usedefault=True)}
    varying_names = []
    for v in varying:
        aname = '{}_map'.format(v)
        desc = 'Path to input {} map'.format(v)
        attrs[aname] = File(exists=True, desc='Path to %s map' % v)
        varying_names.append(aname)
    attrs['varying_names'] = varying_names

    for f in fixed:
        aname = '{}_map'.format(f)
        desc = 'Path to fixed {} map'.format(f)
        args = '--{}=%s'.format(f)
        attrs[aname] = File(argstr=args, exists=True, desc=desc)

    for idx, o in enumerate(out_files):
        aname = '{}_file'.format(o)
        desc = 'Output {} file'.format(o)
        attrs[aname] = File(argstr='%s', mandatory=True,
                            position=idx, desc=desc)
    if extras:
        for k, v in extras.items():
            attrs[k] = v

    T = type(name + 'SimInputSpec', (InputSpec,), attrs)
    return T

################################# Output Specs #################################
#
# QUIT output specifications follow a consistent pattern, so these are defined
# as functions that return a type. Closest thing to C++ templates I could find


def FitOutputSpec(prefix, parameters):
    attrs = {}
    for p in parameters:
        pname = '{}_map'.format(p)
        fname = '{}_{}.nii.gz'.format(prefix, p)
        desc = 'Path to {}'.format(p)
        attrs[pname] = File(fname, desc=desc, usedefault=True)
    attrs['rmse_map'] = File('{}_rmse.nii.gz'.format(
        prefix), desc='Path to RMSE', usedefault=True)
    T = type(prefix + 'FitOutputSpec', (TraitedSpec,), attrs)
    return T


def SimOutputSpec(name, files=None):
    if files is None:
        files = ['out']
    attrs = {}
    for f in files:
        aname = '{}_file'.format(f)
        desc = 'Path to {} image file'.format(f)
        attrs[aname] = File(desc=desc)
    T = type(name + 'SimOutputSpec', (TraitedSpec,), attrs)
    return T

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
            fname = '_' + self.__class__.__name__ + '_input.json'
            with open(fname, 'w') as outfile:
                json.dump(self._json, outfile, indent=2)
            self.inputs.json = fname
        return super()._parse_inputs(skip)


class FitCommand(BaseCommand):
    """
    Support for standard fitting tools.
    """

    def _list_outputs(self):
        outputs = self.output_spec().get()
        prefixed = {}
        for key, trait in outputs.items():
            prefixed[key] = self._gen_fname(trait, prefix=self.inputs.prefix)
        return prefixed

    def __init__(self, sequence={}, **kwargs):
        self._json = deepcopy(sequence)
        super().__init__(**kwargs)


class SimCommand(BaseCommand):
    """
    Base support for QUIT simulator commands.
    """

    def __init__(self, sequence={}, **kwargs):
        self._json = deepcopy(sequence)
        super().__init__(**kwargs)

    def _parse_inputs(self, skip=None):
        if skip is None:
            skip = []
        for v in self.input_spec.varying_names:
            if isdefined(getattr(self.inputs, v)):
                self._json[v] = getattr(self.inputs, v)
                skip.append(v)
        parsed = super()._parse_inputs(skip)
        return parsed

    def _list_outputs(self):
        inputs = self.inputs.get()
        outputs = self.output_spec().get()
        for k, v in outputs.items():
            outputs[k] = path.abspath(inputs[k])
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
