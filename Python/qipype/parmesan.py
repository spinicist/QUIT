#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Implementation of nipype interfaces for QUIT relaxometry utilities
Requires that the QUIT tools are in your your system path
"""

import json
from copy import deepcopy
from os import path, getcwd
from nipype.interfaces.base import traits, isdefined, CommandLine, CommandLineInputSpec, TraitedSpec, DynamicTraitedSpec, File, PackageInfo
from qipype.base import InputSpec, BaseCommand


####
# Redefine these as in fitting for the moment
####

def FitIS(name, fixed=None, in_files=None, extra=None):
    """
    Input Specification for tools in fitting mode
    """
    if fixed is None:
        fixed = []

    if in_files is None:
        in_files = ['in']

    attrs = {'covar': traits.Bool(decsc='Write out parameter covar images',
                                  argstr='--covar'),
             'residuals': traits.Bool(desc='Write out residuals for each data-point',
                                      argstr='--resids'),
             '__module__': __name__}

    for f in fixed:
        aname = '{}_map'.format(f)
        desc = 'Path to fixed {} map'.format(f)
        args = '--{}=%s'.format(f)
        attrs[aname] = File(argstr=args, exists=True, desc=desc)

    for idx, o in enumerate(in_files):
        aname = '{}_file'.format(o)
        desc = 'Output {} file'.format(o)
        attrs[aname] = File(argstr='%s', mandatory=True,
                            position=idx, desc=desc)
    if extra:
        for k, v in extra.items():
            attrs[k] = v
    T = type(name + 'FitIS', (InputSpec,), attrs)
    return T


def SimIS(name, varying, fixed=None, out_files=None, extra=None):
    """
    Input specification for tools in simulation mode
    """
    if fixed is None:
        fixed = []

    if out_files is None:
        out_files = ['out']

    attrs = {'noise': traits.Float(desc='Noise level to add to simulation',
                                   argstr='--simulate=%f', default_value=0.0, usedefault=True),
             '__module__': __name__}
    varying_names = []
    for v in varying:
        aname = '{}_map'.format(v)
        desc = 'Path to input {} map'.format(v)
        attrs[aname] = File(exists=True, desc='Path to {} map'.format(v))
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
    if extra:
        for k, v in extra.items():
            attrs[k] = v
    T = type(name + 'SimIS', (InputSpec,), attrs)
    return T

################################# Output Specs #################################
#
# QUIT output specifications follow a consistent pattern, so these are defined
# as functions that return a type. Closest thing to C++ templates I could find


def FitOS(name, prefix, varying, derived=None):
    attrs = {'__module__': __name__}
    if derived is None:
        derived = []
    for p in varying + derived:
        pname = '{}_map'.format(p)
        fname = '{}_{}.nii.gz'.format(prefix, p)
        desc = 'Path to {}'.format(p)
        attrs[pname] = File(fname, desc=desc, usedefault=True)
    attrs['rmse_map'] = File('{}_rmse.nii.gz'.format(
        prefix), desc='Path to RMSE', usedefault=True)
    T = type(name + 'FitOS', (TraitedSpec,), attrs)
    return T


def SimOS(name, files=None):
    if files is None:
        files = ['out']
    attrs = {'__module__': __name__}
    for f in files:
        aname = '{}_file'.format(f)
        desc = 'Path to {} image file'.format(f)
        attrs[aname] = File(desc=desc)
    T = type(name + 'SimOS', (TraitedSpec,), attrs)
    return T


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
        BaseCommand.__init__(self, **kwargs)


class SimCommand(BaseCommand):
    """
    Base support for QUIT simulator commands.
    """

    def __init__(self, sequence={}, **kwargs):
        self._json = deepcopy(sequence)
        BaseCommand.__init__(self, **kwargs)

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


def Command(toolname, cmd, file_prefix, varying,
            derived=None, fixed=None, files=None, extra=None,
            init=None):
    fit_ispec = FitIS(toolname,
                      fixed=fixed,
                      in_files=files,
                      extra=extra)
    sim_ispec = SimIS(toolname,
                      varying=varying,
                      fixed=fixed,
                      out_files=files,
                      extra=extra)

    fit_ospec = FitOS(toolname, file_prefix,
                      varying=varying, derived=derived)
    sim_ospec = SimOS(toolname, files=files)

    fit_attrs = {'_cmd': cmd, 'input_spec': fit_ispec,
                 'output_spec': fit_ospec}
    sim_attrs = {'_cmd': cmd, 'input_spec': sim_ispec,
                 'output_spec': sim_ospec}
    if init:
        fit_attrs['__init__'] = init
        sim_attrs['__init__'] = init

    fit_cmd = type(toolname, (FitCommand,), fit_attrs)
    sim_cmd = type(toolname + 'Sim', (SimCommand,), sim_attrs)

    return (fit_cmd, sim_cmd, fit_ispec, fit_ospec, sim_ispec, sim_ospec)


Transient, TransientSim, TransientFitIS, TransientFitOS, TransientSimIS, TransientSimOS = Command(
    'Transient', 'qi transient', 'PT',
    varying=['M0', 'T1', 'T2'],
    fixed=['B1'])
TransientB1, TransientB1Sim, TransientB1FitIS, TransientB1FitOS, TransientB1SimIS, TransientB1SimOS = Command(
    'TransientB1', 'qi transient --B1', 'PT_B1',
    varying=['M0', 'T1', 'T2', 'B1'])
SteadyStateT1, SteadyStateT1Sim, SteadyStateT1FitIS, SteadyStateT1FitOS, SteadyStateT1SimIS, SteadyStateT1SimOS = Command(
    'SteadyStateT1', 'qi ss', 'SS',
    varying=['M0', 'T1', 'B1'])
SteadyStateT2, SteadyStateT2Sim, SteadyStateT2FitIS, SteadyStateT2FitOS, SteadyStateT2SimIS, SteadyStateT2SimOS = Command(
    'SteadyStateT2', 'qi ss --T2', 'SS',
    varying=['M0', 'T1', 'T2', 'f0', 'B1'])
SteadyStateMT, SteadyStateMTSim, SteadyStateMTFitIS, SteadyStateMTFitOS, SteadyStateMTSimIS, SteadyStateMTSimOS = Command(
    'SteadyStateMT', 'qi ss --MT', 'SS',
    varying=['M0_f', 'M0_b', 'T1_f', 'T2_f', 'T2_b', 'k', 'f0', 'B1'])
