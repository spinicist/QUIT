/*
 *  MPRAGESequence.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "MPRAGESequence.h"

namespace QI {

Eigen::ArrayXcd MPRAGE::signal(std::shared_ptr<Model> m, const Eigen::VectorXd &par) const {
    return m->MPRAGE(par, FA, TR, ETL, k0, eta, TI, TD);
}

Eigen::ArrayXd MPRAGE::weights(const double f0) const {
    return Eigen::ArrayXd::Ones(size()) * 1.0;
}

Eigen::ArrayXcd MP2RAGE::signal(const double M0, const double T1, const double B1, const double eta) const {
    return One_MP2RAGE(FA, TR, ETL, TD, M0, T1, B1, eta);
}

Eigen::ArrayXcd MP3RAGE::signal(const double M0, const double T1, const double B1, const double eta) const {
    return One_MP3RAGE(FA, TR, ETL, TD, M0, T1, B1, eta);
}

} // End namespace QI
