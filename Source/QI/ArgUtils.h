/*
 *  ArgUtils.h
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "args.hxx"
#include "QI/Types.h"

namespace QI {

template<typename T>
T CheckPos(args::Positional<T> &a) {
    if (a) {
        return a.Get();
    } else {
        std::cerr << a.Name() << " was not specified. Use --help to see usage." << std::endl;
        exit(EXIT_FAILURE);
    }
}

template<typename TRegion = typename QI::VolumeF::RegionType>
TRegion RegionOpt(const std::string &a) {
    std::istringstream iss(a);
    typename TRegion::IndexType start;
    typename TRegion::SizeType size;
    for (int i = 0; i < TRegion::ImageDimension; i++) {
        iss >> start[i];
    }
    for (int i = 0; i < TRegion::ImageDimension; i++) {
        iss >> size[i];
    }
    TRegion r;
    r.SetIndex(start);
    r.SetSize(size);
    return r;
}

}