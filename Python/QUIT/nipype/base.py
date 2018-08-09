# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""The quit module provides classes for interfacing with the `QUIT
<https://github.com/spinicist/QUIT>`_ command line tools.  

"""
from __future__ import (print_function, division, unicode_literals,
                        absolute_import)

from glob import glob
import os

from nipype import logging
from nipype.utils.filemanip import fname_presuffix
from nipype.interfaces.base import (traits, isdefined, CommandLine, CommandLineInputSpec,
                    PackageInfo)
from nipype.external.due import BibTeX



class QUITCommandInputSpec(CommandLineInputSpec):
    """
    Base Input Specification for all QUIT Commands
    """

class FSLCommand(CommandLine):
    """Base support for QUIT commands.

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

   
