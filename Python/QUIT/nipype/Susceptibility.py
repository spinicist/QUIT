#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT Susceptibility.
    
Requires that the QUIT tools are in your your system path

By: Emil Ljungberg and Tobias Wood
"""

from __future__ import (print_function, division, unicode_literals,
                        absolute_import)

from nipype.interfaces.base import CommandLineInputSpec, CommandLine, TraitedSpec, File, traits, isdefined
import json
import os
from nipype.interfaces.fsl import BET
from .base import QUITCommand, QUITCommandInputSpec

############################ qi_unwrap_laplace ############################
# < To be implemented > #

############################ qi_unwrap_path ############################
# < To be implemented > #
