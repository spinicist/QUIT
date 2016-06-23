/*
 *  Kernels.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_KERNELS_H
#define QI_KERNELS_H

#include <memory>
#include <iostream>

#include "QI/Util.h"

namespace QI {

class FilterKernel {
public:
    virtual void print(std::ostream &ostr) const = 0;
    virtual double value(const double &r) const = 0;
};

class TukeyKernel : public FilterKernel {
protected:
    double m_a = 0.75;
    double m_q = 0.25;
public:
    TukeyKernel();
    TukeyKernel(std::istream &istr);
    virtual void print(std::ostream &ostr) const override;
    virtual double value(const double &r) const override;
};

class HammingKernel : public FilterKernel {
protected:
    double m_a = 0.5;
    double m_b = 0.5;
public:
    HammingKernel();
    HammingKernel(std::istream &istr);
    virtual void print(std::ostream &ostr) const override;
    virtual double value(const double &r) const override;
};

class GaussKernel : public FilterKernel {
protected:
    double m_sigma = 0.5;
public:
    GaussKernel();
    GaussKernel(std::istream &istr);
    virtual void print(std::ostream &ostr) const override;
    virtual double value(const double &r) const override;
};

class BlackmanKernel : public FilterKernel {
protected:
    double m_alpha, m_a0, m_a1, m_a2;
public:
    BlackmanKernel();
    BlackmanKernel(std::istream &istr);
    virtual void print(std::ostream &ostr) const override;
    virtual double value(const double &r) const override;
};

std::shared_ptr<FilterKernel> ReadKernel(std::istream &istr);
std::ostream& operator<<(std::ostream &ostr, const FilterKernel &k);

} // End namespace QI

#endif // QI_KERNELS_H