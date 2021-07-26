/*
 *  qi_ase_oef.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Dense>

#include "NumericalIntegration.h"

// #define QI_DEBUG_BUILD

#include "Args.h"

#include "FitFunction.h"
#include "FitScaledAuto.h"
#include "ImageIO.h"
#include "Model.h"
#include "ModelFitFilter.h"
#include "MultiEchoSequence.h"
#include "SimulateModel.h"
#include "Util.h"

using namespace std::literals;

constexpr double kappa      = 0.03;                // Conversion factor
constexpr double gyro_gamma = 2 * M_PI * 42.577e6; // Gyromagnetic Ratio
constexpr double delta_X0   = 0.264e-6;            // Susc diff oxy and fully de-oxy blood

template <typename T> struct fcFunctor {
    const T &dw, tau;
    T        operator()(const T u) const {
        const auto Jterm = ceres::BesselJ0(1.5 * dw * tau * u);
        return (2. + u) * sqrt(1. - u) * (1. - Jterm) / (u * u);
    }
};

struct ASEModel : QI::Model<double, double, 4, 0, 1, 3> {
    using SequenceType = QI::MultiEchoSequence;
    const SequenceType &sequence;
    const double        B0, Hct;

    const std::array<const std::string, NV> varying_names{"S0"s, "dT"s, "R2p"s, "DBV"s};
    const std::array<const std::string, ND> derived_names{"Tc"s, "OEF"s, "dHb"s};

    const VaryingArray start{0.98, 0., 5.0, 0.025};
    const VaryingArray bounds_lo{0.1, -0.1, 0.25, 0.001};
    const VaryingArray bounds_hi{2., 0.1, 50.0, 0.5};

    int input_size(const int /* Unused */) const { return sequence.size(); }

    template <typename Derived>
    auto signal(const Eigen::ArrayBase<Derived> &varying, const FixedArray & /* Unused */) const
        -> QI_ARRAY(typename Derived::Scalar) {
        using T        = typename Derived::Scalar;
        const T &  S0  = varying[0];
        const T &  dT  = varying[1];
        const T &  R2p = varying[2];
        const T &  DBV = varying[3];
        const auto dw  = R2p / DBV;

        Eigen::Integrator<T> integrator(5);
        const auto           quad_rule = Eigen::Integrator<T>::GaussKronrod15;
        T                    abs_error{1.e-9};
        T                    rel_error{1.e-9};
        const auto           aTE = (sequence.TE + dT).abs();
        QI_ARRAY(T) fc(sequence.size());
        for (int i = 0; i < sequence.size(); i++) {
            const auto   tau = aTE(i);
            fcFunctor<T> functor{dw, tau};
            fc[i] = (1. / 3.) * integrator.quadratureAdaptive(
                                    functor, T{0.0}, T{1.0}, abs_error, rel_error, quad_rule);
        }
        QI_ARRAY(T) S = S0 * exp(-DBV * fc);
        return S;
    }

    void derived(const VaryingArray &varying,
                 const FixedArray & /* Unused */,
                 DerivedArray &derived) const {
        const auto &R2p = varying[2];
        const auto &DBV = varying[3];

        const auto dw = R2p / DBV;
        const auto Tc = 1.5 * (1. / dw);
        const auto OEF =
            std::clamp(dw / ((4. * M_PI / 3.) * gyro_gamma * B0 * delta_X0 * Hct), 0., 1.);
        const auto Hb  = Hct / kappa;
        const auto dHb = OEF * Hb;

        derived[0] = Tc;
        derived[1] = OEF;
        derived[2] = dHb;
    }
};
using ASEFit = QI::ScaledAutoDiffFit<ASEModel>;

struct ASEFixDBVModel : QI::Model<double, double, 3, 0, 1, 3> {
    using SequenceType = QI::MultiEchoSequence;
    const SequenceType &sequence;
    const double        B0, Hct, DBV;

    const std::array<const std::string, NV> varying_names{"S0"s, "dT"s, "R2p"s};
    const std::array<const std::string, ND> derived_names{"Tc"s, "OEF"s, "dHb"s};

    const VaryingArray start{0.98, 0., 5.0};
    const VaryingArray bounds_lo{0.1, -0.1, 0.25};
    const VaryingArray bounds_hi{2., 0.1, 50.0};

    int input_size(const int /* Unused */) const { return sequence.size(); }

    template <typename Derived>
    auto signal(const Eigen::ArrayBase<Derived> &varying, const FixedArray & /* Unused */) const
        -> QI_ARRAY(typename Derived::Scalar) {
        using T        = typename Derived::Scalar;
        const T &  S0  = varying[0];
        const T &  dT  = varying[1];
        const T &  R2p = varying[2];
        const auto dw  = R2p / DBV;

        Eigen::Integrator<T> integrator(5);
        const auto           quad_rule = Eigen::Integrator<T>::GaussKronrod15;
        T                    abs_error{1.e-9};
        T                    rel_error{1.e-9};
        const auto           aTE = (sequence.TE + dT).abs();
        QI_ARRAY(T) fc(sequence.size());
        for (int i = 0; i < sequence.size(); i++) {
            const auto   tau = aTE(i);
            fcFunctor<T> functor{dw, tau};
            fc[i] = (1. / 3.) * integrator.quadratureAdaptive(
                                    functor, T{0.0}, T{1.0}, abs_error, rel_error, quad_rule);
        }
        QI_ARRAY(T) S = S0 * exp(-DBV * fc);

        return S;
    }

    void derived(const VaryingArray &varying,
                 const FixedArray & /* Unused */,
                 DerivedArray &derived) const {
        const auto &R2p = varying[2];
        const auto  dw  = R2p / DBV;
        const auto  Tc  = 1.5 * (1. / dw);
        const auto  OEF =
            std::clamp(dw / ((4. * M_PI / 3.) * gyro_gamma * B0 * delta_X0 * Hct), 0., 1.);
        const auto Hb  = Hct / kappa;
        const auto dHb = OEF * Hb;

        derived[0] = Tc;
        derived[1] = OEF;
        derived[2] = dHb;
    }
};
using ASEFixDBVFit = QI::ScaledAutoDiffFit<ASEFixDBVModel>;

/*
 * Main
 */
int ase_oef_main(args::Subparser &parser) {
    args::Positional<std::string> input_path(parser, "ASE_FILE", "Input ASE file");

    QI_COMMON_ARGS;
    args::ValueFlag<double> B0(parser, "B0", "Field-strength (Tesla), default 3", {'B', "B0"}, 3.0);
    args::ValueFlag<double> Hct(parser, "HCT", "Hematocrit (default 0.34)", {'h', "Hct"}, 0.34);
    args::ValueFlag<double> DBV(parser, "DBV", "Fix DBV and only fit R2'", {'d', "DBV"}, 0.0);

    parser.Parse();
    json input    = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    auto sequence = input.at("MultiEcho").get<QI::MultiEchoSequence>();

    if (simulate) {
        if (DBV) {
            ASEFixDBVModel model{{}, sequence, B0.Get(), Hct.Get(), DBV.Get()};
            QI::SimulateModel<ASEFixDBVModel, false>(input,
                                                     model,
                                                     {},
                                                     {input_path.Get()},
                                                     mask.Get(),
                                                     verbose,
                                                     simulate.Get(),
                                                     threads.Get(),
                                                     subregion.Get());
        } else {
            ASEModel model{{}, sequence, B0.Get(), Hct.Get()};
            QI::SimulateModel<ASEModel, false>(input,
                                               model,
                                               {},
                                               {input_path.Get()},
                                               mask.Get(),
                                               verbose,
                                               simulate.Get(),
                                               threads.Get(),
                                               subregion.Get());
        }
    } else {
        auto process = [&](auto fit_func) {
            auto fit_filter = QI::ModelFitFilter<decltype(fit_func)>::New(
                &fit_func, verbose, covar, resids, threads.Get(), subregion.Get());
            fit_filter->ReadInputs({QI::CheckPos(input_path)}, {}, mask.Get());
            fit_filter->Update();
            fit_filter->WriteOutputs(prefix.Get() + "ASE_");
        };
        if (DBV) {
            ASEFixDBVModel model{{}, sequence, B0.Get(), Hct.Get(), DBV.Get()};
            ASEFixDBVFit   fit{model};
            process(fit);
        } else {
            ASEModel model{{}, sequence, Hct.Get(), B0.Get()};
            ASEFit   fit{model};
            process(fit);
        }
    }
    return EXIT_SUCCESS;
}
