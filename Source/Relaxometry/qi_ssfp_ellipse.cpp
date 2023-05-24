/*
 *  qi_ellipse.cpp
 *
 *  Copyright (c) 2016, 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Core>

#include "Args.h"
#include "ImageIO.h"
#include "ModelFitFilter.h"
#include "SSFPSequence.h"
#include "SimulateModel.h"
#include "Util.h"

struct EllipseModel : QI::Model<std::complex<double>, double, 5, 0> {
    QI::SSFPSequence const &sequence;

    std::array<const std::string, NV> const varying_names{"G", "a", "b", "theta_0", "phi_rf"};

    QI_ARRAY(std::complex<double>)
    signal(VaryingArray const &v, FixedArray const & /* Unused */) const {
        const double &G = v[0];
        const double &a = v[1];
        const double &b = v[2];

        if (b > 2. * a / (1. + a * a)) {
            return QI_ARRAY(std::complex<double>)::Zero(sequence.PhaseInc.rows());
        }

        const double &theta0  = v[3];
        const double &psi0    = v[4];
        const auto    theta   = theta0 - sequence.PhaseInc;
        const double  psi     = theta0 / 2.0 + psi0;
        const auto    cos_th  = cos(theta);
        const auto    sin_th  = sin(theta);
        const auto    cos_psi = cos(psi);
        const auto    sin_psi = sin(psi);
        QI_ARRAY(std::complex<double>) result(sequence.PhaseInc.rows());
        result.real() =
            G * (cos_psi - a * (cos_th * cos_psi - sin_th * sin_psi)) / (1.0 - b * cos_th);
        result.imag() =
            G * (sin_psi - a * (cos_th * sin_psi + sin_th * cos_psi)) / (1.0 - b * cos_th);
        return result;
    }

    template <typename Derived>
    auto signal(const Eigen::ArrayBase<Derived> &v) const -> QI_ARRAY(typename Derived::Scalar) {
        using T               = typename Derived::Scalar;
        using ArrayXT         = Eigen::Array<T, Eigen::Dynamic, 1>;
        const T &     G       = v[0];
        const T &     a       = v[1];
        const T &     b       = v[2];
        const T &     theta0  = v[3];
        const T &     psi0    = v[4];
        const ArrayXT theta   = theta0 - sequence.PhaseInc;
        const T       psi     = theta0 / 2.0 + psi0;
        const ArrayXT cos_th  = cos(theta);
        const ArrayXT sin_th  = sin(theta);
        const T       cos_psi = cos(psi);
        const T       sin_psi = sin(psi);
        const ArrayXT re_m =
            (cos_psi - a * (cos_th * cos_psi - sin_th * sin_psi)) * G / (1.0 - b * cos_th);
        const ArrayXT im_m =
            (sin_psi - a * (cos_th * sin_psi + sin_th * cos_psi)) * G / (1.0 - b * cos_th);
        ArrayXT result(re_m.rows() + im_m.rows());
        result << re_m, im_m;
        return result;
    }
};

struct EllipseCost {
    EllipseModel const &model;
    QI_ARRAY(std::complex<double>) const data;

    template <typename T> bool operator()(const T *const vin, T *rin) const {
        const Eigen::Map<const QI_ARRAYN(T, EllipseModel::NV)> v(vin);

        // Check if ellipse has gone horizontal
        const T &a = v[1];
        const T &b = v[2];
        if (b < 2. * a / (1. + a * a)) {
            QI_ARRAY(T) m = model.signal(v);
            Eigen::Map<QI_ARRAY(T)> r(rin, data.size() * 2);
            r.head(data.size()) = m.head(data.size()) - data.real();
            r.tail(data.size()) = m.tail(data.size()) - data.imag();
            return true;
        } else {
            /* ELLIPSE COST NOT EVALUATED - ELLIPSE WAS HORIZONTAL */
            return false;
        }
    }
};

struct EllipseFit {
    static const bool Blocked = true;
    static const bool Indexed = false;
    using InputType           = std::complex<double>;
    using OutputType          = double;
    using RMSErrorType        = double;
    using FlagType            = int;
    using ModelType           = EllipseModel;
    ModelType model;

    int input_size(const int /* Unused */) const { return model.sequence.size(); }
    int n_outputs() const { return model.NV; }

    QI::FitReturnType fit(const std::vector<Eigen::ArrayXcd> &inputs,
                          EllipseModel::FixedArray const &    fixed,
                          EllipseModel::VaryingArray &        p,
                          EllipseModel::CovarArray *          cov,
                          double &                            rmse,
                          std::vector<Eigen::ArrayXcd> &      residuals,
                          FlagType &                          iterations,
                          const int /* Unused */) const {
        const double               scale  = inputs[0].abs().maxCoeff();
        const Eigen::ArrayXcd      data   = inputs[0] / scale;
        const std::complex<double> c_mean = data.mean();

        using AutoCost = ceres::AutoDiffCostFunction<EllipseCost, ceres::DYNAMIC, EllipseModel::NV>;
        auto *auto_cost           = new AutoCost(new EllipseCost{model, data}, data.rows() * 2);
        ceres::LossFunction *loss = new ceres::HuberLoss(1.0);
        ceres::Problem       problem;
        problem.AddResidualBlock(auto_cost, loss, p.data());
        const double not_zero = 1.0e-6;
        const double not_one  = 1.0 - not_zero;
        const double max_a    = exp(-model.sequence.TR / 5.0); // Set a sensible maximum on T2
        problem.SetParameterLowerBound(p.data(), 0, not_zero);
        problem.SetParameterUpperBound(p.data(), 0, not_one);
        problem.SetParameterLowerBound(p.data(), 1, not_zero);
        problem.SetParameterUpperBound(p.data(), 1, max_a);
        problem.SetParameterLowerBound(p.data(), 2, not_zero);
        problem.SetParameterUpperBound(p.data(), 2, not_one);
        problem.SetParameterLowerBound(p.data(), 3, -2. * M_PI);
        problem.SetParameterUpperBound(p.data(), 3, 2. * M_PI);
        problem.SetParameterLowerBound(p.data(), 4, -2. * M_PI);
        problem.SetParameterUpperBound(p.data(), 4, 2. * M_PI);
        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
        options.max_num_iterations  = 50;
        options.function_tolerance  = 1e-5;
        options.gradient_tolerance  = 1e-6;
        options.parameter_tolerance = 1e-3;
        options.logging_type        = ceres::SILENT;
        double th0 = 0.0, psi0 = 0.0, best_cost = std::numeric_limits<double>::infinity();
        for (const auto &th0_try : {-M_PI, 0., M_PI}) {
            const double psi0_try = arg(c_mean / std::polar(1.0, th0_try / 2));
            p << abs(c_mean), 0.5, 0.5, th0_try, psi0_try;
            double cost = 0.0;
            problem.Evaluate(ceres::Problem::EvaluateOptions(), &cost, NULL, NULL, NULL);
            if (cost < best_cost) {
                best_cost = cost;
                th0       = th0_try;
                psi0      = psi0_try;
            }
        }
        p << abs(c_mean), 0.5, 0.5, th0, psi0;
        ceres::Solve(options, &problem, &summary);
        if (!summary.IsSolutionUsable()) {
            return {false, summary.FullReport()};
        }
        iterations = summary.iterations.size();

        Eigen::ArrayXcd const rs  = (data - model.signal(p, fixed));
        double const          var = rs.abs().square().sum();
        rmse                      = sqrt(var / data.rows()) * scale;
        if (residuals.size() > 0) {
            residuals[0] = rs * scale;
        }
        if (cov) {
            QI::GetModelCovariance<ModelType>(problem, p, var / (data.rows() - ModelType::NV), cov);
        }

        p[0] *= scale;
        p[3] = std::fmod(p[3] + 3 * M_PI, 2 * M_PI) - M_PI;
        p[4] = std::fmod(p[4] + 3 * M_PI, 2 * M_PI) - M_PI;
        return {true, ""};
    }
};

int ssfp_ellipse_main(args::Subparser &parser) {
    args::Positional<std::string> sequence_path(parser, "sequence_FILE", "Input sequence file");
    QI_COMMON_ARGS;
    args::ValueFlag<char> algorithm(
        parser, "ALGO", "Choose algorithm (h)yper/(d)irect, default d", {'a', "algo"}, 'd');
    parser.Parse();
    QI::CheckPos(sequence_path);
    QI::Log(verbose, "Reading sequence information");
    json input    = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    auto sequence = input.at("SSFP").get<QI::SSFPSequence>();
    QI::Log(verbose, "{}", sequence);
    EllipseModel model{{}, sequence};
    if (simulate) {
        QI::SimulateModel<EllipseModel, false>(input,
                                               model,
                                               {},
                                               {sequence_path.Get()},
                                               mask.Get(),
                                               verbose,
                                               simulate.Get(),
                                               threads.Get(),
                                               subregion.Get());
    } else {
        EllipseFit fit{model};
        auto       fit_filter =
            QI::ModelFitFilter<EllipseFit>::New(
                &fit, verbose, covar, resids, threads.Get(), subregion.Get());
        fit_filter->ReadInputs({sequence_path.Get()}, {}, mask.Get());
        fit_filter->SetBlocks(fit_filter->GetInput(0)->GetNumberOfComponentsPerPixel() /
                              sequence.size());
        fit_filter->Update();
        fit_filter->WriteOutputs(prefix.Get() + "ES_");
        QI::Log(verbose, "Finished.");
    }
    return EXIT_SUCCESS;
}
