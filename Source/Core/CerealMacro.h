/*
 *  CerealUtil.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_CEREALMACRO_H
#define QI_CEREALMACRO_H

#include <type_traits>
#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include "Macro.h"

namespace QI {

template<typename Archive, typename T>
inline void cereal_load(Archive &ar, T &x, const std::string &name, const std::string &location) {
    try {
        ar(cereal::make_nvp(name, x));
    } catch (cereal::RapidJSONException &e) {
        std::cerr << "Error loading parameter " << name << " during: " << location << std::endl;
        exit(EXIT_FAILURE);
    };
}

template<typename Archive, typename T>
inline void cereal_load_degrees(Archive &ar, T &x, const std::string &name, const std::string &location) {
    try {
        T temp;
        ar(cereal::make_nvp(name, temp));
        x = temp * M_PI / 180.;
    } catch (cereal::RapidJSONException &e) {
        std::cerr << "Error loading parameter " << name << " during: " << location << std::endl;
        exit(EXIT_FAILURE);
    };
}

template<typename Archive, typename T>
inline void cereal_save(Archive &ar, const T &x, const std::string &name, const std::string &location) {
    try {
        ar(cereal::make_nvp(name, x));
    } catch (cereal::RapidJSONException &e) {
        std::cerr << "Error saving parameter " << name << " during: " << location << std::endl;
        exit(EXIT_FAILURE);
    };
}

template<typename Archive, typename T>
inline void cereal_save_degrees(Archive &ar, const T &x, const std::string &name, const std::string &location) {
    try {
        T temp = x * 180. / M_PI;
        ar(cereal::make_nvp(name, temp));

    } catch (cereal::RapidJSONException &e) {
        std::cerr << "Error saving parameter " << name << " during: " << location << std::endl;
        exit(EXIT_FAILURE);
    };
}

} // End namespace QI

#define QI_CLOAD( AR, X ) QI::cereal_load( AR, X, #X, QI_LOCATION );
#define QI_CLOAD_NAME( AR, X, NAME ) QI::cereal_load( AR, X, NAME, QI_LOCATION );
#define QI_CLOAD_DEGREES( AR, X ) QI::cereal_load_degrees( AR, X, #X, QI_LOCATION );

#define QI_CSAVE( AR, X ) QI::cereal_save( AR, X, #X, QI_LOCATION );
#define QI_CSAVE_NAME( AR, X, NAME ) QI::cereal_save( AR, X, NAME, QI_LOCATION );
#define QI_CSAVE_DEGREES( AR, X ) QI::cereal_save_degrees( AR, X, #X, QI_LOCATION );

#endif // QI_CEREALMACRO_H
