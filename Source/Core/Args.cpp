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

namespace QI {

void ParseArgs(args::ArgumentParser &parser, int argc, char **argv, const args::Flag &verbose) {
    try {
        parser.ParseCLI(argc, argv);
        QI::Info(verbose, "Starting {} {}", argv[0], QI::GetVersion());
    } catch (args::Help) {
        fmt::print("{}\n", parser);
        exit(EXIT_SUCCESS);
    } catch (args::ParseError e) {
        QI::Fail("{}\n{}", parser, e.what());
    } catch (args::ValidationError e) {
        QI::Fail("{}\n{}", parser, e.what());
    }
}

void ParseArgs(args::ArgumentParser &parser,
               int                   argc,
               char **               argv,
               const args::Flag &    verbose,
               args::ValueFlag<int> &threads) {
    try {
        parser.ParseCLI(argc, argv);
        QI::Info(verbose, "Starting {} {}", argv[0], QI::GetVersion());
        QI::Log(verbose, "Max threads = {}", threads.Get());
        itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads(threads.Get());
    } catch (args::Help) {
        fmt::print("{}\n", parser);
        exit(EXIT_SUCCESS);
    } catch (args::ParseError e) {
        QI::Fail("{}\n{}", parser, e.what());
    } catch (args::ValidationError e) {
        QI::Fail("{}\n{}", parser, e.what());
    }
}

} // End namespace QI
