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

args::Group    global_group("GLOBAL OPTIONS");
args::HelpFlag help(global_group, "HELP", "Show this help message", {'h', "help"});
args::Flag     verbose(global_group, "VERBOSE", "Talk more", {'v', "verbose"});
