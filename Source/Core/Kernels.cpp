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

#include "Kernels.h"

namespace QI {
TukeyKernel::TukeyKernel() {}
TukeyKernel::TukeyKernel(std::istream &istr) {
    if (!istr.eof()) {
        std::string nextValue;
        std::getline(istr, nextValue, ',');
        m_a = stod(nextValue);
        std::getline(istr, nextValue, ',');
        m_q = stod(nextValue);
    }
}
void TukeyKernel::print(std::ostream &ostr) const {
    ostr << "Tukey," << m_a << "," << m_q;
}
double TukeyKernel::value(const Eigen::Array3d &pos, const Eigen::Array3d &sz, const Eigen::Array3d &sp) const {
    const double r = sqrt(((pos / sz).square()).sum() / 3);
    const double v = (r <= (1 - m_a)) ? 1 : 0.5*((1+m_q)+(1-m_q)*cos(M_PI*(r - (1 - m_a))/m_a));
    return v;
}

HammingKernel::HammingKernel() {}
HammingKernel::HammingKernel(std::istream &istr) {
    if (!istr.eof()) {
        std::string nextValue;
        std::getline(istr, nextValue, ',');
        m_a = stod(nextValue);
        std::getline(istr, nextValue, ',');
        m_b = stod(nextValue);
    }
}
void HammingKernel::print(std::ostream &ostr) const {
    ostr << "Hamming," << m_a << "," << m_b;
}
double HammingKernel::value(const Eigen::Array3d &pos, const Eigen::Array3d &sz, const Eigen::Array3d &sp) const {
    const double r = sqrt(((pos / sz).square()).sum() / 3);
    const double v = m_a - m_b*cos(M_PI*(1.+r));
    return v;
}

GaussKernel::GaussKernel() {}
GaussKernel::GaussKernel(std::istream &istr) {
    if (!istr.eof()) {
        std::string nextValue;
        std::getline(istr, nextValue, ',');
        if (istr) { // Still more values
            m_fwhm[0] = stod(nextValue);
            std::getline(istr, nextValue, ',');
            m_fwhm[1] = stod(nextValue);
            std::getline(istr, nextValue);
            m_fwhm[2] = stod(nextValue);
        } else {
            m_fwhm = Eigen::Array3d::Ones() * stod(nextValue);
        }
    }
}
void GaussKernel::print(std::ostream &ostr) const {
    ostr << "Gauss," << m_fwhm.transpose();
}

double GaussKernel::value(const Eigen::Array3d &pos, const Eigen::Array3d &sz, const Eigen::Array3d &sp) const {
    static const double M = 2. * sqrt(2.*log(2.)) / M_PI;
    const Eigen::Array3d sigma_k = M * sz * sp / m_fwhm;
    const double r2 = (pos/sigma_k).square().sum();
    const double v = exp(-r2/2.);
    return v;
}

BlackmanKernel::BlackmanKernel() {
    calc_constants();
}
BlackmanKernel::BlackmanKernel(std::istream &istr) {
    if (!istr.eof()) {
        std::string nextValue;
        std::getline(istr, nextValue, ',');
        m_alpha = stod(nextValue);
    }
    calc_constants();
}
void BlackmanKernel::calc_constants() {
    m_a0 = (1. - m_alpha) / 2.;
    m_a1 = 1. / 2.;
    m_a2 = m_alpha / 2.;
}
void BlackmanKernel::print(std::ostream &ostr) const {
    ostr << "Blackman," << m_alpha;
}
double BlackmanKernel::value(const Eigen::Array3d &pos, const Eigen::Array3d &sz, const Eigen::Array3d &sp) const {
    const double r = sqrt(((pos / sz).square() / 3).sum());
    const double v = m_a0 - m_a1*cos(M_PI*(1.+r)) + m_a2*cos(2.*M_PI*(1.+r));
    return v;
}

RectKernel::RectKernel() {}
RectKernel::RectKernel(std::istream &istr) {
    if (!istr.eof()) {
        std::string nextValue;
        std::getline(istr, nextValue, ','); m_dim = stoi(nextValue);
        std::getline(istr, nextValue, ','); m_width = stoi(nextValue);
        std::getline(istr, nextValue, ','); m_val_inside = stod(nextValue);
        std::getline(istr, nextValue); m_val_outside = stod(nextValue);
    }
    if (m_dim > 2) {
        QI_EXCEPTION("Dimension for Rectangular filter must be less than 3");
    }
}
void RectKernel::print(std::ostream &ostr) const {
    ostr << "FixFSE," << m_dim << "," << m_width << "," << m_val_inside << "," << m_val_outside;
}
double RectKernel::value(const Eigen::Array3d &pos, const Eigen::Array3d &sz, const Eigen::Array3d &sp) const {
    if (fabs(pos[m_dim]) > m_width) {
        return m_val_outside;
    } else {
        return m_val_inside;
    }
}

FixFSEKernel::FixFSEKernel() {
}
FixFSEKernel::FixFSEKernel(std::istream &istr) {
    if (!istr.eof()) {
        std::string nextValue;
        std::getline(istr, nextValue, ','); m_dim = stoi(nextValue);
        std::getline(istr, nextValue, ','); m_etl = stoi(nextValue);
        std::getline(istr, nextValue, ','); m_kzero = stoi(nextValue);
        std::getline(istr, nextValue, ','); m_te1 = stod(nextValue);
        std::getline(istr, nextValue, ','); m_esp = stod(nextValue);
        std::getline(istr, nextValue); m_T2 = stod(nextValue);
    }
    if (m_dim > 2) {
        QI_EXCEPTION("Dimension for FixFSE filter must be less than 3");
    }
}
void FixFSEKernel::print(std::ostream &ostr) const {
    ostr << "FixFSE," << m_dim << "," << m_etl << "," << m_kzero
         << "," << m_te1 << "," << m_esp << "," << m_T2;
}
double FixFSEKernel::value(const Eigen::Array3d &pos, const Eigen::Array3d &sz, const Eigen::Array3d &sp) const {
    const int dim = abs(m_dim);
    const int dir = m_dim > 0 ? 1 : -1;
    const int n_trains = 2 * sz[dim] / m_etl;
    const int n_echo = floor((dir*pos[dim] + sz[dim] - (m_etl / 2) + 2) / n_trains);
    // At this point, center of kspace is at m_etl / 2. Shift to make it kzero
    const int n_shifted = n_echo - (m_etl / 2) + m_kzero;
    // Wrap negative echoes to end of train
    const int n_final = (n_shifted < 0) ? n_shifted + m_etl : n_shifted;
    const double E2 = exp((m_te1 + n_final*m_esp)/m_T2);
    return E2;
}

std::shared_ptr<FilterKernel> ReadKernel(const std::string &str) {
    std::shared_ptr<FilterKernel> newKernel = nullptr;
    std::istringstream iss(str);
    std::string filterName;
    std::getline(iss, filterName, ',');
    if (filterName == "Tukey") {
        newKernel = std::make_shared<TukeyKernel>(iss);
    } else if (filterName == "Hamming") {
        newKernel = std::make_shared<HammingKernel>(iss);
    } else if (filterName == "Gauss") {
        newKernel = std::make_shared<GaussKernel>(iss);
    } else if (filterName == "Blackman") {
        newKernel = std::make_shared<BlackmanKernel>(iss);
    } else if (filterName == "Rectangle") {
        newKernel = std::make_shared<RectKernel>(iss);
    } else if (filterName == "FixFSE") {
        newKernel = std::make_shared<FixFSEKernel>(iss);
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