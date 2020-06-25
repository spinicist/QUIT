/*
 *  Args.h
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#pragma once

#include "ImageTypes.h"
#include "Log.h"
#include "args.hxx"

namespace QI {

template <typename T> T CheckPos(args::Positional<T> &a) {
    if (a) {
        return a.Get();
    } else {
        QI::Fail("{} was not specified. Use --help to see usage.", a.Name());
    }
}

template <typename T> std::vector<T> CheckList(args::PositionalList<T> &a) {
    if (a) {
        return a.Get();
    } else {
        QI::Fail("No values of {} specified. Use --help to see usage.", a.Name());
    }
}

template <typename T> T CheckValue(args::ValueFlag<T> &v) {
    if (v) {
        return v.Get();
    } else {
        QI::Fail("{} was not specified but is required. Use --help to see usage.", v.Name());
    }
}

template <typename TArray, const int size> void ArrayArg(const std::string &a, TArray &array) {
    std::istringstream iss(a);
    std::string        el;
    for (int i = 0; i < size; i++) {
        std::getline(iss, el, ',');
        if (!iss) {
            QI::Fail("Failed to read from array argument: {}", a);
        }
        array[i] = std::stoi(el);
    }
}

template <typename TArray, const int size> void ArrayArgF(const std::string &a, TArray &array) {
    std::istringstream iss(a);
    std::string        el;
    for (int i = 0; i < size; i++) {
        std::getline(iss, el, ',');
        if (!iss) {
            QI::Fail("Failed to read from array argument: {}", a);
        }
        array[i] = std::stof(el);
    }
}

} // End namespace QI

extern args::Group    global_group;
extern args::HelpFlag help;
extern args::Flag     verbose;

#define QI_COMMON_ARGS                                                                         \
    args::Flag resids(parser, "RESIDS", "Write point residuals", {'r', "resids"});             \
    args::Flag covar(                                                                          \
        parser, "COVAR", "Write out covariance matrix (CoV and Corr) images", {"covar"});      \
                                                                                               \
    args::ValueFlag<int>   threads(parser,                                                     \
                                 "THREADS",                                                  \
                                 "Use N threads (default=hardware limit or $QUIT_THREADS)",  \
                                 {'T', "threads"},                                           \
                                 QI::GetDefaultThreads());                                   \
    args::ValueFlag<float> simulate(                                                           \
        parser, "SIMULATE", "Simulate sequence (argument is noise level)", {"simulate"}, 0.0); \
    args::ValueFlag<std::string> mask(                                                         \
        parser, "MASK", "Only process voxels within the given mask", {'m', "mask"});           \
    args::ValueFlag<std::string> subregion(                                                    \
        parser,                                                                                \
        "SUBREGION",                                                                           \
        "Process voxels in a block from I,J,K with size SI,SJ,SK",                             \
        {'s', "subregion"});                                                                   \
    args::ValueFlag<std::string> prefix(                                                       \
        parser, "PREFIX", "Add a prefix to output filenames", {'o', "out"});                   \
    args::ValueFlag<std::string> json_file(                                                    \
        parser, "JSON", "Read JSON from file instead of stdin", {"json"});
