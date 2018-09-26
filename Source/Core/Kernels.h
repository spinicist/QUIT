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

#include "Eigen/Core"
#include "Util.h"

namespace QI
{

class FilterKernel {
protected:

public:
    virtual void print(std::ostream &ostr) const = 0;
    virtual double value(const Eigen::Array3d &pos, const Eigen::Array3d &sz, const Eigen::Array3d &sp) const = 0;
    virtual ~FilterKernel() = default;
};

class TukeyKernel : public FilterKernel
{
  protected:
    double m_a = 0.75;
    double m_q = 0.25;

  public:
    TukeyKernel();
    TukeyKernel(std::istream &istr);
    virtual void print(std::ostream &ostr) const override;
    virtual double value(const Eigen::Array3d &pos, const Eigen::Array3d &sz, const Eigen::Array3d &sp) const override;
};

class HammingKernel : public FilterKernel
{
  protected:
    double m_a = 0.5;
    double m_b = 0.5;

  public:
    HammingKernel();
    HammingKernel(std::istream &istr);
    virtual void print(std::ostream &ostr) const override;
    virtual double value(const Eigen::Array3d &pos, const Eigen::Array3d &sz, const Eigen::Array3d &sp) const override;
};

class GaussKernel : public FilterKernel
{
  protected:
    Eigen::Array3d m_fwhm;

  public:
    GaussKernel();
    GaussKernel(std::istream &istr);
    virtual void print(std::ostream &ostr) const override;
    virtual double value(const Eigen::Array3d &pos, const Eigen::Array3d &sz, const Eigen::Array3d &sp) const override;
};

class BlackmanKernel : public FilterKernel
{
  protected:
    double m_alpha = 0.16, m_a0, m_a1, m_a2;
    void calc_constants(); // Helper function for constructors

  public:
    BlackmanKernel();
    BlackmanKernel(std::istream &istr);
    virtual void print(std::ostream &ostr) const override;
    virtual double value(const Eigen::Array3d &pos, const Eigen::Array3d &sz, const Eigen::Array3d &sp) const override;
};

class RectKernel : public FilterKernel {
protected:
    int m_dim = 0, m_width = 8;
    double m_val_inside = 1, m_val_outside = 0;
  
public:
    RectKernel();
    RectKernel(std::istream &istr);
    virtual void print(std::ostream &ostr) const override;
    virtual double value(const Eigen::Array3d &pos,
                         const Eigen::Array3d &sz,
                         const Eigen::Array3d &sp) const override;
};

class FixFSEKernel : public FilterKernel {
protected:
    int m_dim = 0, m_etl = 16, m_kzero = 6;
    double m_te1 = 0, m_esp = 8, m_T2 = 40;
  
public:
    FixFSEKernel();
    FixFSEKernel(std::istream &istr);
    virtual void print(std::ostream &ostr) const override;
    virtual double value(const Eigen::Array3d &pos,
                         const Eigen::Array3d &sz,
                         const Eigen::Array3d &sp) const override;
};

std::shared_ptr<FilterKernel> ReadKernel(const std::string &str);
std::ostream &operator<<(std::ostream &ostr, const FilterKernel &k);

} // End namespace QI

#endif // QI_KERNELS_H