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

#ifndef QI_ARGS_H
#define QI_ARGS_H

#include "args.hxx"
#include "Macro.h"
#include "ImageTypes.h"
#include "Util.h"

namespace QI {

void ParseArgs(args::ArgumentParser &parser, int argc, char **argv, const args::Flag &verbose) {
    try {
        parser.ParseCLI(argc, argv);
        if (verbose) std::cout << "Starting " << argv[0] << " " << QI::GetVersion() << std::endl;
    } catch (args::Help) {
        std::cout << parser;
        exit(EXIT_SUCCESS);
    } catch (args::ParseError e) {
        QI_FAIL(e.what() << std::endl << parser);
    } catch (args::ValidationError e) {
        QI_FAIL(e.what() << std::endl << parser);
    }
}

template<typename T>
T CheckPos(args::Positional<T> &a) {
    if (a) {
        return a.Get();
    } else {
        QI_FAIL(a.Name() << " was not specified. Use --help to see usage.");
    }
}

template<typename T>
std::vector<T> CheckList(args::PositionalList<T> &a) {
    if (a) {
        return a.Get();
    } else {
        QI_FAIL("No values of " << a.Name() << " specified. Use --help to see usage.");
    }
}

template<typename T>
T CheckValue(args::ValueFlag<T> &v) {
    if (v) {
        return v.Get();
    } else {
        QI_FAIL(v.Name() << " was not specified. Use --help to see usage.");
    }
}

template<typename TArray, const int size>
void ArrayArg(const std::string &a, TArray &array) {
    std::istringstream iss(a);
    std::string el;
    for (int i = 0; i < size; i++) {
        std::getline(iss, el, ',');
        if (!iss) {
            QI_FAIL("Failed to read array argument from string: " << a);
        }
        array[i] = std::stoi(el);
    }
}

template<typename TRegion = typename QI::VolumeF::RegionType>
TRegion RegionArg(const std::string &a) {
    std::istringstream iss(a);
    std::string el;
    typename TRegion::IndexType start;
    typename TRegion::SizeType size;
    for (int i = 0; i < TRegion::ImageDimension; i++) {
        std::getline(iss, el, ',');
        start[i] = std::stoi(el);
    }
    for (int i = 0; i < TRegion::ImageDimension; i++) {
        std::getline(iss, el, ',');
        size[i] = std::stoi(el);
    }
    TRegion r;
    r.SetIndex(start);
    r.SetSize(size);
    return r;
}

} // End namespace QI

#endif // QI_ARGS_H