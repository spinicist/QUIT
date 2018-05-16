#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Helper functions for QUIT nipype interfaces
"""

import os
import json

def parse_param_dict(name, spec, value):
    if name == 'param_dict':
        with open('_tmp_input.json', 'w') as outfile:
            json.dump(value, outfile)
        return "< _tmp_input.json"

    return super()._format_arg(name, spec, value)