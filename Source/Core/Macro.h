/*
 *  Macro.h
 *  Part of the QUantitative Image Toolbox
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

/*
 * The following macros are adapted from ITK for use in QUIT files which don't include any of ITK
 */

#ifndef QI_MACRO_H
#define QI_MACRO_H

#define FMT_DEPRECATED_OSTREAM

#include "fmt/format.h"
#include "fmt/ostream.h"
#include <Eigen/Core>

#if defined(_WIN32) && !defined(__MINGW32__)
#define QI_LOCATION __FUNCSIG__
#elif defined(__GNUC__)
#define QI_LOCATION __PRETTY_FUNCTION__
#else
#define QI_LOCATION __FUNCTION__
#endif

#define QI_ARRAY(T) Eigen::Array<T, Eigen::Dynamic, 1>
#define QI_ARRAYN(T, N) Eigen::Array<T, N, 1>

#ifdef QI_DEBUG_BUILD
#define QI_DBMSG(x) fmt::print(x);
#define QI_DB(x) fmt::print("{}: {}\n", #x, x);
#define QI_DBMAT(x) fmt::print("{}:\n{}\n", #x, x);
#define QI_DBVEC(x) fmt::print("{}: {}\n", #x, x.transpose());
#define QI_DBSTL(x) fmt::print("{}: {}\n", #x, fmt::join(x, ", "));
#else
#define QI_DBMSG(x)
#define QI_DB(x)
#define QI_DBMAT(x)
#define QI_DBVEC(x)
#define QI_DBVECT(x)
#define QI_DBSTL(x)
#endif

#define QI_ForwardNewMacro(x)                                      \
    template <class... Args> static Pointer New(Args &&... args) { \
        Pointer smartPtr = ::itk::ObjectFactory<x>::Create();      \
        if (smartPtr == nullptr) {                                 \
            smartPtr = new x(std::forward<Args>(args)...);         \
        }                                                          \
        smartPtr->UnRegister();                                    \
        return smartPtr;                                           \
    }

#endif // QI_MACRO_H