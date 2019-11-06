/*
 *  Util.h
 *  Part of the QUantitative Image Toolbox
 *
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QUIT_UTIL_H
#define QUIT_UTIL_H

#include <Eigen/Core>
#include <functional>
#include <map>
#include <mutex>
#include <random>
#include <string>
#include <vector>

#include "itkVector.h"

#include "Log.h"

/*
 * Print a std::vector
 */
template <typename T> inline std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec) {
    os << "(";
    if (vec.size() > 0) {
        auto it = vec.begin();
        std::cout << *it;
        for (it++; it != vec.end(); it++) {
            std::cout << ", " << *it;
        }
    }
    os << ")";
    return os;
}

namespace QI {

int GetDefaultThreads(); //!< Return the number of threads in the $QUIT_THREADS environment variable
const std::string &GetVersion(); //!< Return the version of the QI library
const std::string &OutExt();     //!< Return the extension stored in $QUIT_EXT

std::string StripExt(const std::string &filename); //!< Remove the extension from a filename
std::string GetExt(const std::string &filename);   //!< Return the extension from a filename with .
std::string Basename(const std::string &path);     //!< Return only the filename part of a path
std::vector<size_t> SortedUniqueIndices(Eigen::ArrayXd const &x); //!< For splines
std::vector<int>    IntsFromString(const std::string &s); // !!< Ints from comma-separated string
std::mt19937_64::result_type RandomSeed();                //!< Thread-safe random seed

/*
 * Helper function to calculate the volume of a voxel in an image
 */
template <typename TImg> double VoxelVolume(const itk::SmartPointer<TImg> &img) {
    double vox_volume = img->GetSpacing()[0];
    if (TImg::ImageDimension < 3) {
        QI::Fail("Tried to find the volume of a <3D image");
    }
    for (int v = 1; v < 3; v++)
        vox_volume *= img->GetSpacing()[v];
    return vox_volume;
}

/*
 * Helper function to clamp between two values
 */
template <typename T> T Clamp(const T &value, const T &low, const T &high) {
    if (value > low) {
        if (value < high) {
            return value;
        } else {
            return high;
        }
    } else {
        return low;
    }
}

template <typename T, unsigned int D>
Eigen::Array<T, D, 1> Eigenify(const itk::Vector<T, D> &itk_vector) {
    Eigen::Array<T, D, 1> eigen_vector;
    for (size_t i = 0; i < D; i++) {
        eigen_vector[i] = itk_vector[i];
    }
    return eigen_vector;
}

template <typename TRegion> TRegion RegionFromString(const std::string &a) {
    auto                        ints = IntsFromString(a);
    typename TRegion::IndexType start;
    typename TRegion::SizeType  size;
    for (size_t i = 0; i < TRegion::ImageDimension; i++) {
        start[i] = ints[i];
        size[i]  = ints[i + TRegion::ImageDimension];
    }
    TRegion r;
    r.SetIndex(start);
    r.SetSize(size);
    return r;
}

} // namespace QI

#endif // QUIT_UTIL_H
