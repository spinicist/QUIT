/*
 *  Kernels.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "QI/Kernels.h"

namespace QI {
TukeyKernel::TukeyKernel() {}
TukeyKernel::TukeyKernel(std::istream &istr) {
    std::string nextValue;
    std::getline(istr, nextValue, ',');
    m_a = stod(nextValue);
    std::getline(istr, nextValue, ',');
    m_q = stod(nextValue);
}
void TukeyKernel::print(std::ostream &ostr) const {
    ostr << "Tukey," << m_a << "," << m_q << std::endl;
}
double TukeyKernel::value(const double &r) const {
    const double v = (r <= (1 - m_a)) ? 1 : 0.5*((1+m_q)+(1-m_q)*cos((M_PI/m_a)*(r - 1 + m_a)));
    return v;
}

HammingKernel::HammingKernel() {}
HammingKernel::HammingKernel(std::istream &istr) {
    std::string nextValue;
    std::getline(istr, nextValue, ',');
    m_a = stod(nextValue);
    std::getline(istr, nextValue, ',');
    m_b = stod(nextValue);
}
void HammingKernel::print(std::ostream &ostr) const {
    ostr << "Hamming," << m_a << "," << m_b << std::endl;
}
double HammingKernel::value(const double &r) const {
    const double v = m_a - m_b*cos(M_PI*(1.+r));
    return v;
}

GaussKernel::GaussKernel() {}
GaussKernel::GaussKernel(std::istream &istr) {
    std::string nextValue;
    std::getline(istr, nextValue, ',');
    m_sigma = stod(nextValue);
}
void GaussKernel::print(std::ostream &ostr) const {
    ostr << "Gauss," << m_sigma << std::endl;
}
double GaussKernel::value(const double &r) const {
    const double v = exp(-pow(r/m_sigma,2)/2.);
    return v;
}

BlackmanKernel::BlackmanKernel() {
    m_alpha = 0.16;
    m_a0 = (1. - m_alpha) / 2.;
    m_a1 = 1. / 2.;
    m_a2 = m_alpha / 2.;
}
BlackmanKernel::BlackmanKernel(std::istream &istr) {
    std::string nextValue;
    std::getline(istr, nextValue, ',');
    m_alpha = stod(nextValue);
    m_a0 = (1. - m_alpha) / 2.;
    m_a1 = 1. / 2.;
    m_a2 = m_alpha / 2.;
}
void BlackmanKernel::print(std::ostream &ostr) const {
    ostr << "Blackman," << m_alpha << std::endl;
}
double BlackmanKernel::value(const double &r) const {
    const double v = m_a0 - m_a1*cos(M_PI*(1.+r)) + m_a2*cos(2.*M_PI*(1.+r));
    return v;
}

std::shared_ptr<FilterKernel> ReadKernel(std::istream &istr) {
    std::shared_ptr<FilterKernel> newKernel = nullptr;
    std::string filterName;
    std::getline(istr, filterName, ',');
    if (filterName == "Tukey") {
        newKernel = std::make_shared<TukeyKernel>(istr);
    } else if (filterName == "Hamming") {
        newKernel = std::make_shared<HammingKernel>(istr);
    } else if (filterName == "Gauss") {
        newKernel = std::make_shared<GaussKernel>(istr);
    } else if (filterName == "Blackman") {
        newKernel = std::make_shared<BlackmanKernel>(istr);
    } else {
        QI_EXCEPTION("Unknown filter type");
    }
    return newKernel;
}

std::ostream& operator<<(std::ostream &ostr, const FilterKernel &k) {
    k.print(ostr);
    return ostr;
}

} // End namespace QI