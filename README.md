![Logo](Docs/logo.png)

[![Build Status](https://dev.azure.com/spinicist/QUIT/_apis/build/status/QUIT-CI?branchName=master)](https://dev.azure.com/spinicist/QUIT/_build/latest?definitionId=4&branchName=master)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00656/status.svg)](https://doi.org/10.21105/joss.00656)
[![DOI](https://zenodo.org/badge/37066948.svg)](https://zenodo.org/badge/latestdoi/37066948)

Credit / Blame / Contact - Tobias Wood - tobias.wood@kcl.ac.uk

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
 
If you find the tools useful the author would love to hear from you.

This is the updated version of QUIT based on ITK http://www.itk.org
The previous version is here http://github.com/spinicist/old_quit

## Brief Description

A collection of programs for processing quantitative MRI data, originally DESPOT
but now a wider variety of techniques.

## Thanks

Many thanks to Samuel Hurley for many, many suggestions.
Thanks to Anna Combes, Amy McDowell and Sjoerd Vos for finding bugs in early
versions.

## Documentation

Documentation exists in Markdown format in the Docs folder, and is also
available at http://spinicist.github.io/QUIT/

## Installation

Pre-compiled binaries are provided for Linux and Mac OS X in a .tar.gz archive
from http://github.com/spinicist/QUIT/releases.

Download the correct archive for your platform, untar it and then ensure that
the binaries can be found via your `PATH` environment variable. The Linux
binaries are built with Ubuntu 14.04 with GCC 6. If you need to run on an older
version of Linux with a previous version of `glibc` then you will need to
compile from source.

For instructions on how to compile from source, please see the developer page
in the documentation.

## Usage

Example scripts for mcDESPOT processing are provided in the installation
archive. These print usage instructions if you call them with no arguments.
They take a set of filenames as input, and you will need to modify the scripts
with your particular flip-angles and TRs.

Each product has some basic usage instructions that will be printed with either
the -h or --help options, e.g. `qimcdespot -h`. The majority of programs take
the input filenames as arguments on the command-line, and in addition will read
an input text file from `stdin`. For further details, see the [documentation](https://spinicist.github.io/QUIT), which is also available in the `Docs` folder.

## Getting Help

If you can't find an answer to a problem in the documentation or help strings, you can open an [issue](https://github.com/spinicist/QUIT/issues), post a question on [Neurostars](https://neurostars.org) or find the main developer on Twitter ([@spinicist](https://twitter.com/spinicist)).
