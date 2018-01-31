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

#include <string>
#include <map>
#include <vector>
#include <random>
#include <functional>
#include <mutex>
#include <vector>

#include <Eigen/Core>

#include "itkCommand.h"

#include "Macro.h"

/*
 * Print a std::vector
 */
template<typename T>
inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
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

const std::string &GetVersion();                    //!< Return the version of the QI library
const std::string &OutExt();                        //!< Return the extension stored in $QUIT_EXT
std::string StripExt(const std::string &filename);  //!< Remove the extension from a filename
std::string GetExt(const std::string &filename);    //!< Return the extension from a filename (including .)
std::string Basename(const std::string &path);      //!< Return only the filename part of a path
std::mt19937_64::result_type RandomSeed();          //!< Thread-safe random seed
unsigned long long Choose(unsigned long long n, unsigned long long k); //!< From Knuth, surprised this isn't in STL
Eigen::ArrayXXd ReadArrayFile(const std::string &path);

class GenericMonitor : public itk::Command {
public:
    typedef GenericMonitor          Self;
    typedef itk::Command            Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    itkTypeMacro(GenericMonitor, Superclass);
    itkNewMacro(Self);
protected:
    GenericMonitor() {}
public:
    void Execute(itk::Object *caller, const itk::EventObject &event) ITK_OVERRIDE;
    void Execute(const itk::Object *object, const itk::EventObject &event) ITK_OVERRIDE;
};

/*
 * Helper function to calculate the volume of a voxel in an image
 */
template<typename TImg>
double VoxelVolume(const itk::SmartPointer<TImg> &img) {
    double vox_volume = img->GetSpacing()[0];
    if (TImg::ImageDimension < 3) {
        QI_EXCEPTION("Tried to find the volume of a <3D image");
    }
    for (int v = 1; v < 3; v++)
        vox_volume *= img->GetSpacing()[v];
    return vox_volume;
}

/*
 * Helper function to clamp between two values
 */
template<typename T> T Clamp(const T &value, const T &low, const T &high) {
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

} // End namespace QUIT

#endif // QUIT_UTIL_H
