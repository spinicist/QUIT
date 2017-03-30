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

#include <sstream>

#if defined( _WIN32 ) && !defined( __MINGW32__ )
    #define QI_LOCATION __FUNCSIG__
#elif defined( __GNUC__ )
    #define QI_LOCATION __PRETTY_FUNCTION__   
#else
    #define QI_LOCATION __FUNCTION__
#endif

#define QI_EXCEPTION( x )                               \
{                                                       \
    std::ostringstream message;                         \
    message << QI_LOCATION << std::endl                 \
            << __FILE__ << ":" << __LINE__ << std::endl \
            << x << std::endl;                          \
    throw(std::runtime_error(message.str()));                 \
}

#define QI_FAIL( x )             \
{                                \
    std::cerr << x << std::endl; \
    exit(EXIT_FAILURE);          \
}

#endif // QI_MACRO_H