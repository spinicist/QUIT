/*
 *  Polynomial.h
 *  Part of the QUantitative Image Toolbox
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_POLYNOMIAL_H
#define QI_POLYNOMIAL_H

#include <Eigen/Core>

namespace QI {

// From Knuth, surprised this isn't in STL
unsigned long long Choose(unsigned long long n, unsigned long long k) {
    if (k > n)
        return 0;

    unsigned long long r = 1;
    for (unsigned long long d = 1; d <= k; ++d) {
        r *= n--;
        r /= d;
    }
    return r;
}

template <int D> class Polynomial {
  protected:
    static const int Dimension = D;
    int              m_order;
    Eigen::ArrayXd   m_coeffs;

  public:
    Polynomial() : m_order(0), m_coeffs(1) { m_coeffs.setConstant(0); }

    Polynomial(const int o) : m_order(o), m_coeffs(QI::Choose(o + Dimension, o)) {
        m_coeffs.setConstant(0);
    }

    Polynomial(const Polynomial &p) : m_order(p.m_order), m_coeffs(p.m_coeffs) {}

    int                   order() const { return m_order; }
    const Eigen::ArrayXd &coeffs() const { return m_coeffs; }
    void                  setCoeffs(const Eigen::ArrayXd &c) { m_coeffs = c; }

    int            nterms() const { return m_coeffs.rows(); }
    Eigen::ArrayXd terms(const Eigen::Vector3d &p) {
        Eigen::ArrayXd                          ts(m_coeffs.rows());
        int                                     it = 0;
        Eigen::Matrix<double, Dimension + 1, 1> all;
        all << 1, p;
        std::function<void(double, int, int)> orderLoop = [&](double t, int o, int start) -> void {
            if (o == m_order) {
                ts[it++] = t;
            } else {
                for (int i = start; i < Dimension + 1; i++) {
                    orderLoop(t * all[i], o + 1, i);
                }
            }
        };
        orderLoop(1, 0, 0);
        return ts;
    }

    Eigen::VectorXd values(const Eigen::Vector3d &p) { return terms(p) * m_coeffs; }

    double value(const Eigen::Vector3d &p) { return values(p).sum(); }

    std::string get_terms() const {
        std::vector<std::string> terms(m_coeffs.rows(), "a");
        std::ostringstream       output;
        std::string              list = "1";
        char                     dim  = 'a';
        for (int i = 0; i < Dimension; i++) {
            list.append(std::string(1, dim));
            dim++;
        }
        int it = 0;

        std::function<void(std::string, int, int)> orderLoop =
            [&](std::string term, int o, int start) -> void {
            if (o == m_order) {
                terms[it++] = term;
            } else {
                for (int i = start; i < Dimension + 1; i++) {
                    orderLoop(term + list[i], o + 1, i);
                }
            }
        };
        orderLoop(std::string(""), 0, 0);
        output << terms[0];
        for (int i = 1; i < m_coeffs.rows(); i++) {
            output << " + " << terms[i];
        }
        return output.str();
    }
};

} // End namespace QI

#endif // QI_POLYNOMIAL_H