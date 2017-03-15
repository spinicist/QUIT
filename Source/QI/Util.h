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

#include "QI/Macro.h"
#include "QI/Args.h"
#include "QI/Types.h"

namespace QI {

const std::string &GetVersion();                    //!< Return the version of the QI library
const std::string &OutExt();                        //!< Return the extension stored in $QUIT_EXT
std::string StripExt(const std::string &filename);  //!< Remove the extension from a filename
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

} // End namespace QUIT

#endif // QUIT_UTIL_H
