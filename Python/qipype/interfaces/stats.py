#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT Stats.
    
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

############################ qi_glmsetup ############################
# < To be implemented > #

############################ qi_glmcontrasts ############################
# < To be implemented > #

############################ qi_rois ############################
# < To be implemented > #