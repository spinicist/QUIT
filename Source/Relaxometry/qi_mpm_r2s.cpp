/*
 *  qi_mpm_r2s.cpp
 *
 *  Copyright (c) 2019 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Core>

#include "Args.h"
#include "ImageIO.h"
#include "Model.h"
#include "ModelFitFilter.h"
#include "MultiEchoSequence.h"
#include "SimulateModel.h"
#include "Util.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

struct MPMModel : QI::Model<double, double, 4, 0, 3> {
    QI::MultiEchoSequence &pdw_s, &t1w_s, &mtw_s;

    double noise = 0.;

    std::array<std::string, NV> const varying_names{"R2s", "S0_PDw", "S0_T1w", "S0_MTw"};

    VaryingArray const lo{1e-6, 1e-6, 1e-6, 1e-6};
    VaryingArray const hi{1e4, 1e2, 1e2, 1e2}; // Signal values will be scaled

    int output_size(int const o) {
        switch (o) {
        case 0:
            return pdw_s.size();
        case 1:
            return t1w_s.size();
        case 2:
            return mtw_s.size();
        default:
            QI::Fail("Incorrect output requested {}", o);
        }
    }
    size_t num_outputs() const { return 3; }

    template <typename Derived>
    auto pdw_signal(const Eigen::ArrayBase<Derived> &v) const
        -> QI_ARRAY(typename Derived::Scalar) {
        using T     = typename Derived::Scalar;
        T const &R2 = v[0];
        T const &PD = v[1]; // S_PDw
        return PD * exp(-pdw_s.TE * R2);
    }

    template <typename Derived>
    auto t1w_signal(const Eigen::ArrayBase<Derived> &v) const
        -> QI_ARRAY(typename Derived::Scalar) {
        using T     = typename Derived::Scalar;
        T const &R2 = v[0];
        T const &PD = v[2]; // S_T1w
        return PD * exp(-pdw_s.TE * R2);
    }

    template <typename Derived>
    auto mtw_signal(const Eigen::ArrayBase<Derived> &v) const
        -> QI_ARRAY(typename Derived::Scalar) {
        using T     = typename Derived::Scalar;
        T const &R2 = v[0];
        T const &PD = v[3]; // S_MTw
        return PD * exp(-pdw_s.TE * R2);
    }

    auto signals(const QI_ARRAYN(double, NV) & v, const QI_ARRAYN(double, NF) & /* Unused */) const
        -> std::vector<QI_ARRAY(double)> {
        return {pdw_signal(v), t1w_signal(v), mtw_signal(v)};
    }
};

struct PDwCost {
    MPMModel const &model;
    QI_ARRAY(double) const data;

    template <typename T> bool operator()(const T *const vin, T *rin) const {
        Eigen::Map<QI_ARRAYN(T, MPMModel::NV) const> const v(vin);

        Eigen::Map<QI_ARRAY(T)> r(rin, data.rows());

        if (model.noise > 0.) {
            r = (data.square() - model.noise) - model.pdw_signal(v).square();
        } else {
            r = data - model.pdw_signal(v);
        }
        return true;
    }
};

struct T1wCost {
    MPMModel const &model;
    QI_ARRAY(double) const data;

    template <typename T> bool operator()(const T *const vin, T *rin) const {
        Eigen::Map<QI_ARRAYN(T, MPMModel::NV) const> const v(vin);

        Eigen::Map<QI_ARRAY(T)> r(rin, data.rows());

        if (model.noise > 0.) {
            r = (data.square() - model.noise) - model.t1w_signal(v).square();
        } else {
            r = data - model.t1w_signal(v);
        }
        return true;
    }
};

struct MTwCost {
    MPMModel const &model;
    QI_ARRAY(double) const data;

    template <typename T> bool operator()(const T *const vin, T *rin) const {
        Eigen::Map<QI_ARRAYN(T, MPMModel::NV) const> const v(vin);

        Eigen::Map<QI_ARRAY(T)> r(rin, data.rows());

        if (model.noise > 0.) {
            r = (data.square() - model.noise) - model.mtw_signal(v).square();
        } else {
            r = data - model.mtw_signal(v);
        }
        return true;
    }
};

struct MPMFit {
    static const bool Blocked = false;
    static const bool Indexed = false;
    using InputType           = double;
    using OutputType          = double;
    using RMSErrorType        = double;
    using FlagType            = int;
    using ModelType           = MPMModel;
    ModelType model;

    int input_size(const int i) const {
        switch (i) {
        case 0:
            return model.pdw_s.size();
        case 1:
            return model.t1w_s.size();
        case 2:
            return model.mtw_s.size();
        default:
            QI::Fail("Invalid input size = {}", i);
        }
    }
    int n_outputs() const { return model.NV; }

    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs,
                          const Eigen::ArrayXd & /* Unused */,
                          ModelType::VaryingArray &    v,
                          ModelType::CovarArray *      cov,
                          RMSErrorType &               rmse,
                          std::vector<Eigen::ArrayXd> &residuals,
                          FlagType &                   iterations) const {
        double scale = std::max({inputs[0].maxCoeff(), inputs[1].maxCoeff(), inputs[2].maxCoeff()});
        if (scale < std::numeric_limits<double>::epsilon()) {
            v    = ModelType::VaryingArray::Zero();
            rmse = 0.0;
            return {false, "Maximum data value was zero or less"};
        }
        Eigen::ArrayXd const pdw_data = inputs[0] / scale;
        Eigen::ArrayXd const t1w_data = inputs[1] / scale;
        Eigen::ArrayXd const mtw_data = inputs[2] / scale;
        v << 20., 1., 1., 1.; // R2s, S_PDw, S_T1w, S_MTw
        ceres::Problem problem;
        using AutoPDwType = ceres::AutoDiffCostFunction<PDwCost, ceres::DYNAMIC, ModelType::NV>;
        using AutoT1wType = ceres::AutoDiffCostFunction<T1wCost, ceres::DYNAMIC, ModelType::NV>;
        using AutoMTwType = ceres::AutoDiffCostFunction<MTwCost, ceres::DYNAMIC, ModelType::NV>;
        auto *pdw_cost    = new AutoPDwType(new PDwCost{model, pdw_data}, model.pdw_s.size());
        auto *t1w_cost    = new AutoT1wType(new T1wCost{model, t1w_data}, model.t1w_s.size());
        auto *mtw_cost    = new AutoMTwType(new MTwCost{model, mtw_data}, model.mtw_s.size());
        ceres::LossFunction *loss = new ceres::HuberLoss(1.0);
        problem.AddResidualBlock(pdw_cost, loss, v.data());
        problem.AddResidualBlock(t1w_cost, loss, v.data());
        problem.AddResidualBlock(mtw_cost, loss, v.data());
        for (int i = 0; i < ModelType::NV; i++) {
            problem.SetParameterLowerBound(v.data(), i, model.lo[i]);
            problem.SetParameterUpperBound(v.data(), i, model.hi[i]);
        }
        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
        options.max_num_iterations  = 50;
        options.function_tolerance  = 1e-5;
        options.gradient_tolerance  = 1e-6;
        options.parameter_tolerance = 1e-4;
        options.logging_type        = ceres::SILENT;
        ceres::Solve(options, &problem, &summary);
        if (!summary.IsSolutionUsable()) {
            return {false, summary.FullReport()};
        }
        iterations = summary.iterations.size();

        Eigen::ArrayXd const pdw_resid = pdw_data - model.pdw_signal(v);
        Eigen::ArrayXd const t1w_resid = t1w_data - model.t1w_signal(v);
        Eigen::ArrayXd const mtw_resid = mtw_data - model.mtw_signal(v);
        if (residuals.size() > 0) {
            residuals[0] = pdw_resid * scale;
            residuals[1] = t1w_resid * scale;
            residuals[2] = mtw_resid * scale;
        }
        double const var =
            pdw_resid.square().sum() + t1w_resid.square().sum() + mtw_resid.square().sum();
        int const dsize = model.pdw_s.size() + model.t1w_s.size() + model.mtw_s.size();
        if (cov) {
            QI::GetModelCovariance<MPMModel>(problem, v, var / (dsize - ModelType::NV), cov);
        }
        rmse      = sqrt(var / dsize);
        v.tail(3) = v.tail(3) * scale; // Multiply signals/proton densities back up
        return {true, ""};
    }
};

/*
 * Main
 */
int mpm_r2s_main(args::Subparser &parser) {
    args::Positional<std::string> pdw_path(parser, "PDw", "Input multi-echo PD-weighted file");
    args::Positional<std::string> t1w_path(parser, "T1w", "Input multi-echo T1-weighted file");
    args::Positional<std::string> mtw_path(parser, "MTw", "Input multi-echo MT-weighted file");
    args::ValueFlag<double>       rician_noise(
        parser, "RICIAN", "Mean squared noise level for Rician correction", {"rician"}, 0.);
    QI_COMMON_ARGS;
    parser.Parse();
    QI::CheckPos(pdw_path);
    QI::CheckPos(t1w_path);
    QI::CheckPos(mtw_path);

    QI::Log(verbose, "Reading sequence parameters");
    json input = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    QI::MultiEchoSequence pdw_seq(input["PDw"]), t1w_seq(input["T1w"]), mtw_seq(input["MTw"]);

    MPMModel model{{}, pdw_seq, t1w_seq, mtw_seq, rician_noise.Get()};
    MPMFit   mpm_fit{model};
    if (simulate) {
        QI::SimulateModel<MPMModel, true>(input,
                                          model,
                                          {},
                                          {pdw_path.Get(), t1w_path.Get(), mtw_path.Get()},
                                          mask.Get(),
                                          verbose,
                                          simulate.Get(),
                                          subregion.Get());
    } else {
        auto fit_filter =
            QI::ModelFitFilter<MPMFit>::New(&mpm_fit, verbose, covar, resids, subregion.Get());
        fit_filter->ReadInputs({pdw_path.Get(), t1w_path.Get(), mtw_path.Get()}, {}, mask.Get());
        fit_filter->Update();
        fit_filter->WriteOutputs(prefix.Get() + "MPM_");
    }
    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
