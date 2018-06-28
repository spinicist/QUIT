/*
 *  EigenCeral.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_EIGENCEREAL_H
#define QI_EIGENCEREAL_H

#include <Eigen/Dense>
#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>

#include "Macro.h"

namespace Eigen {
    template<typename Archive, typename T, cereal::traits::EnableIf<cereal::traits::is_text_archive<Archive>::value> = cereal::traits::sfinae>
    inline void load(Archive &ar, Array<T, Dynamic, 1> &v) {
        cereal::size_type n_rows;
        ar(cereal::make_size_tag(n_rows));
        if (static_cast<cereal::size_type>(v.rows()) != n_rows) {
            v.resize(n_rows, 1);
        }

        for(cereal::size_type i = 0; i < n_rows; i++) {
            ar(v[i]);
        }
    }

    template<typename Archive, typename T, cereal::traits::EnableIf<cereal::traits::is_text_archive<Archive>::value> = cereal::traits::sfinae>
    inline void save(Archive &ar, const Array<T, Dynamic, 1> &v) {
        ar(cereal::make_size_tag(v.rows()));
        for(auto i = 0; i < v.rows(); ++i) {
            ar(v[i]);
        }
    }

    template<typename Archive, typename T, int D,
             cereal::traits::EnableIf<cereal::traits::is_text_archive<Archive>::value> = cereal::traits::sfinae>
    inline void load(Archive &ar, Array<T, D, 1> &v) {
        cereal::size_type sz = D;
        ar(cereal::make_size_tag(sz));
        for (size_t i = 0; i < D; i++) {
            ar(v[i]);
        }
    }

    template<typename Archive, typename T, int D,
             cereal::traits::EnableIf<cereal::traits::is_text_archive<Archive>::value> = cereal::traits::sfinae>
    inline void save(Archive &ar, const Array<T, D, 1> &v) {
        cereal::size_type sz = D;
        ar(cereal::make_size_tag(sz));
        for (size_t i = 0; i < D; i++) {
            ar(v[i]);
        }
    }

} // end namespace Eigen

#endif // End QI_EIGENCEREAL_H