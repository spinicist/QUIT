/*
 *  Lineshape.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef LINESHAPE_H
#define LINESHAPE_H

#include "JSON.h"
#include "Macro.h"
#include "NumericalIntegration.h"
#include "ceres/cubic_interpolation.h"
#include <Eigen/Core>
#include <functional>
#include <memory>
#include <string>

namespace QI {

enum class Lineshapes { Gaussian, Lorentzian, SuperLorentzian, Interpolated };

template <typename T> QI_ARRAY(T) Gaussian(const Eigen::ArrayXd &f0, const T T2b) {
    return sqrt(1.0 / (2.0 * M_PI)) * T2b * exp(-pow(2.0 * M_PI * f0 * T2b, 2.0) / 2.0);
}

template <typename T> QI_ARRAY(T) Lorentzian(const Eigen::ArrayXd &f0, const T T2b) {
    return T2b / (1.0 + pow(2.0 * M_PI * f0 * T2b, 2.0));
}

template <typename T> struct SLFunctor {
    const T &     T2b;
    const double &f0;
    T             operator()(const T u) const {
        const auto uterm = abs(T2b / (3.0 * u * u - 1.0));
        return sqrt(2.0 / M_PI) * uterm * exp(-2.0 * pow(2.0 * M_PI * f0 * uterm, 2.0));
    }
};

template <typename T> QI_ARRAY(T) SuperLorentzian(const Eigen::ArrayXd &df0, const T T2b) {
    Eigen::Integrator<T> integrator(200);
    const auto           quad_rule = Eigen::Integrator<T>::GaussKronrod61;
    T                    abs_error{0.0};
    T                    rel_error{Eigen::NumTraits<double>::epsilon() * 50.0};
    QI_ARRAY(T) vals(df0.rows());
    for (auto i = 0; i < df0.rows(); i++) {
        SLFunctor<T> sl_functor{T2b, df0[i]};
        vals[i] = integrator.quadratureAdaptive(
            sl_functor, T{0.0}, T{1.0}, abs_error, rel_error, quad_rule);
    }
    return vals;
}

struct InterpLineshape {
    double                                                           T2_nominal = 1e-6;
    double                                                           freq_min, freq_step;
    int                                                              freq_count;
    Eigen::ArrayXd                                                   values;
    std::shared_ptr<ceres::Grid1D<double>>                           grid;
    std::shared_ptr<ceres::CubicInterpolator<ceres::Grid1D<double>>> interpolator;

    InterpLineshape(const double          freq_min,
                    const double          freq_step,
                    const int             freq_count,
                    const Eigen::ArrayXd &vals,
                    const double          T2b);
    InterpLineshape(double const T2b, Eigen::ArrayXd &freqs, Eigen::ArrayXd &vals);

    template <typename T> QI_ARRAY(T) operator()(const Eigen::ArrayXd &f, const T T2) const {
        const auto scale = T2 / T2_nominal;
        const auto sf    = (f.abs() * scale - freq_min) / freq_step;
        QI_ARRAY(T) interp_vals(f.rows());
        for (auto i = 0; i < f.rows(); i++) {
            if (sf[i] < 0.0) {
                interp_vals[i] = values[0] * scale;
            } else if (sf[i] > (freq_count - 1.0)) {
                interp_vals[i] = values[freq_count - 1] * scale;
            } else {
                interpolator->Evaluate(sf[i], &interp_vals[i]);
                interp_vals[i] *= scale;
            }
        }
        return interp_vals;
    }

    template <typename T> T operator()(double const &f, const T T2) const {
        const auto scale = T2 / T2_nominal;
        const auto sf    = (std::abs(f) * scale - freq_min) / freq_step;

        if (sf < 0.0) {
            return values[0] * scale;
        } else if (sf > (freq_count - 1.0)) {
            return values[freq_count - 1] * scale;
        } else {
            T val;
            interpolator->Evaluate(sf, &val);
            val *= scale;
            return val;
        }
    }
};

} // End namespace QI

namespace nlohmann {
template <> struct adl_serializer<QI::InterpLineshape> {
    static QI::InterpLineshape from_json(const json &j);
    static void                to_json(json &j, QI::InterpLineshape t);
};
} // namespace nlohmann

#endif // LINESHAPE_H
