# -*- coding: utf-8 -*-
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
