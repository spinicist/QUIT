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
protected:
    int m_hx = 0, m_hy = 0, m_hz = 0;
    double radius(const int x, const int y, const int z) const;
public:
    void setSize(const int sx, const int sy, const int sz);
    virtual void print(std::ostream &ostr) const = 0;
    virtual double value(const int x, const int y, const int z) const = 0;
};

class TukeyKernel : public FilterKernel {
protected:
    double m_a = 0.75;
    double m_q = 0.25;
public:
    TukeyKernel();
    TukeyKernel(std::istream &istr);
    virtual void print(std::ostream &ostr) const override;
    virtual double value(const int x, const int y, const int z) const override;
};

class HammingKernel : public FilterKernel {
protected:
    double m_a = 0.5;
    double m_b = 0.5;
public:
    HammingKernel();
    HammingKernel(std::istream &istr);
    virtual void print(std::ostream &ostr) const override;
    virtual double value(const int x, const int y, const int z) const override;
};

class GaussKernel : public FilterKernel {
protected:
    double m_a = 0.5;
public:
    GaussKernel();
    GaussKernel(std::istream &istr);
    virtual void print(std::ostream &ostr) const override;
    virtual double value(const int x, const int y, const int z) const override;
};

std::shared_ptr<FilterKernel> ReadKernel(std::istream &istr);
std::ostream& operator<<(std::ostream &ostr, const FilterKernel &k);

} // End namespace QI

#endif // QI_KERNELS_H