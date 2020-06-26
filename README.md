![Logo](Docs/logo.png)

[![Build](https://github.com/spinicist/QUIT/workflows/Build/badge.svg)](https://github.com/spinicist/QUIT/actions)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00656/status.svg)](https://doi.org/10.21105/joss.00656)
[![DOI](https://zenodo.org/badge/37066948.svg)](https://zenodo.org/badge/latestdoi/37066948)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/spinicist/qipype-binder/master?filepath=T1_mapping.ipynb)


Credit / Blame / Contact - Tobias Wood - tobias.wood@kcl.ac.uk

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
 
If you find the tools useful the author would love to hear from you.

## Brief Description

A collection of programs for processing quantitative MRI data, originally DESPOT
but now a wider variety of techniques.

## Thanks

Many thanks to Samuel Hurley, Erika Raven, Anna Combes, Amy McDowell and 
Sjoerd Vos.

## Documentation

Documentation is available at http://quit.readthedocs.io (and the in Docs folder
in .rst format).

## Installation

### Executable

Pre-compiled executables are provided for Linux and Mac OS X in a .tar.gz 
archive from http://github.com/spinicist/QUIT/releases. Download the archive and 
extract it with `tar -xzf qi-platform.tar.gz`. Then, move the resulting `qi` 
file to somewhere on your `$PATH`, for instance `/usr/local/bin`. That's it.

- MacOS Catalina or higher users should use `curl` to download the binary, i.e. 
  type `curl -L https://github.com/spinicist/QUIT/releases/download/v3.0/qi-macos.tar.gz`
  This is because Safari now sets the quarantine attribute of all downloads,
  which prevents them being run as the binary is unsigned. It is possible to 
  remove the quarantine flag with `xattr`, but downloading with `curl` is more 
  straightforward.
- The Linux executable is compiled on Ubuntu 16.04 with GLIBC version 2.3 and a 
  statically linked libc++. This means it will hopefully run on most modern 
  Linux distributions. Let me know if it doesn't.

For instructions on how to compile from source, please see the developer page
in the documentation.

### Python

QUIT now comes with `nipype` wrappers. Install with `pip install qipype`.

## Usage

QUIT comes as a single executable file with multiple commands, similar to `git` 
or `bart`. Type `qi` to see a list of all the available commands. The majority 
of commands require an input 4D Nifti file containing the image data to fit, 
and a `.json` file containing the sequence parameters (e.g. `TR`). See 
examples for each command at https://quit.readthedocs.io.

An example iPython notebook is available here in a separate repo https://github.com/spinicist/qipype-binder. Workflows for CEST and MPM are provided in `qipype.workflows`.

## Getting Help

If you can't find an answer to a problem in the documentation or help strings, 
you can open an [issue](https://github.com/spinicist/QUIT/issues), or find the 
developer on Twitter ([@spinicist](https://twitter.com/spinicist)).
