/*
 *  Args.cpp
 *
 *  Copyright (c) 2019 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "Args.h"
#include "Util.h"

#include "itkMultiThreaderBase.h"

args::Group          global_group("GLOBAL OPTIONS");
args::HelpFlag       help(global_group, "HELP", "Show this help message", {'h', "help"});
args::Flag           verbose(global_group, "VERBOSE", "Talk more", {'v', "verbose"});
args::ValueFlag<int> threads(global_group,
                             "THREADS",
                             "Use N threads (default=hardware limit or $QUIT_THREADS)",
                             {'T', "threads"},
                             QI::GetDefaultThreads());

void Parse(args::Subparser &parser) {
    parser.Parse();
    QI::SetVerbose(verbose.Get());
    itk::MultiThreaderBase::SetGlobalDefaultNumberOfThreads(threads.Get());
}