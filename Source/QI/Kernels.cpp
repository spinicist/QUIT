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

double FilterKernel::radius(const int x, const int y, const int z) const {
    const double rx = fmod(static_cast<double>(x)/m_hx + 1.0, 2.0) - 1.0;
    const double ry = fmod(static_cast<double>(y)/m_hy + 1.0, 2.0) - 1.0;
    const double rz = fmod(static_cast<double>(z)/m_hz + 1.0, 2.0) - 1.0;
    const double r = sqrt((rx*rx + ry*ry + rz*rz) / 3);
    return r;
}

void FilterKernel::setSize(const int sx, const int sy, const int sz) {
    m_hx = sx/2; m_hy = sy/2; m_hz = sz/2;
}

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
double TukeyKernel::value(const int x, const int y, const int z) const {
    const double r = radius(x, y, z);
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
double HammingKernel::value(const int x, const int y, const int z) const {
    const double r = radius(x, y, z);
    const double v = m_a + m_b*cos(2.*M_PI*r);
    return v;
}

GaussKernel::GaussKernel() {}
GaussKernel::GaussKernel(std::istream &istr) {
    std::string nextValue;
    std::getline(istr, nextValue, ',');
    m_a = stod(nextValue);
}
void GaussKernel::print(std::ostream &ostr) const {
    ostr << "Gauss," << m_a << std::endl;
}
double GaussKernel::value(const int x, const int y, const int z) const {
    const double r = radius(x, y, z);
    const double v = exp(-m_a*r*r);
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