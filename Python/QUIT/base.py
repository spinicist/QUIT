# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""The quit module provides classes for interfacing with the `QUIT
<https://github.com/spinicist/QUIT>`_ command line tools.  

Heavily inspired by the FSL base class
"""

from __future__ import (print_function, division, unicode_literals,
                        absolute_import)

from glob import glob
import os
import json

from nipype import logging
from nipype.utils.filemanip import fname_presuffix
from nipype.interfaces.base import (traits, isdefined, CommandLine, CommandLineInputSpec,
                                    PackageInfo)
from nipype.external.due import BibTeX


class QUITCommandInputSpec(CommandLineInputSpec):
    """
    Base Input Specification for all QUIT Commands
    """

    # Inputs that are common to all program
    verbose = traits.Bool(desc='Print more information', argstr='-v')
    environ = {'QUIT_EXT': 'NIFTI_GZ'}


class QUITCommand(CommandLine):
    """
    Base support for QUIT commands.
    """

    input_spec = QUITCommandInputSpec
    _output_type = None

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

    def _add_prefix(self, full_path):
        """
        Add prefix to file in full_path
        """
        if self.inputs.prefix:
            p, f = os.path.split(full_path)
            return os.path.join(p, self.inputs.prefix + f)
        else:
            return full_path

    def _gen_fname(self,
                   basename,
                   cwd=None,
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
            cwd = os.getcwd()
        if suffix is None:
            suffix = ''
        fname = fname_presuffix(
            basename, suffix=suffix, use_ext=True, newpath=cwd)
        return fname

    def _process_params(self, name, spec, value):
        """
        Make parameter dictionary into a .json file for input to interface
        """
        if name == 'param_dict':
            with open('_tmp_input.json', 'w') as outfile:
                json.dump(value, outfile)
            return "< _tmp_input.json"

        return super()._format_arg(name, spec, value)


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
