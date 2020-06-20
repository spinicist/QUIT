/*
 *  despot1hifi.cpp
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright (c) 2011-2013 Tobias Wood.
 *
 *  Based on code by Sean Deoni
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "ceres/ceres.h"
#include <Eigen/Core>
#include <array>

#include "Args.h"
#include "FitFunction.h"
#include "ImageIO.h"
#include "Model.h"
#include "ModelFitFilter.h"
#include "OnePoolSignals.h"
#include "SequenceGroup.h"
#include "SimulateModel.h"
#include "Util.h"

using namespace std::literals;

struct HIFIModel : QI::Model<double, double, 3, 0, 2> {
    QI::SPGRSequence   spgr;
    QI::MPRAGESequence mprage;

    const std::array<const std::string, 3> varying_names{"PD"s, "T1"s, "B1"s};

    VaryingArray bounds_lo = QI_ARRAYN(double, NV)::Constant(1.0e-4);
    VaryingArray bounds_hi =
        QI_ARRAYN(double, NV)::Constant(std::numeric_limits<double>::infinity());

    size_t num_outputs() const { return 2; }
    int    output_size(int i) {
        if (i == 0) {
            return spgr.size();
        } else if (i == 1) {
            return mprage.size();
        } else {
            QI::Fail("Invalid output size: {}", i);
        }
    }

    template <typename Derived>
    auto spgr_signal(const Eigen::ArrayBase<Derived> &v) const
        -> QI_ARRAY(typename Derived::Scalar) {
        return QI::SPGRSignal(v[0], v[1], v[2], spgr);
    }

    template <typename Derived>
    auto mprage_signal(const Eigen::ArrayBase<Derived> &v) const
        -> QI_ARRAY(typename Derived::Scalar) {
        return QI::MPRAGESignal(v[0], v[1], v[2], mprage);
    }

    auto signals(const QI_ARRAYN(double, NV) & v, const QI_ARRAYN(double, NF) & /* Unused */) const
        -> std::vector<QI_ARRAY(double)> {
        return {spgr_signal(v), mprage_signal(v)};
    }
};

struct HIFISPGRCost {
    const HIFIModel &model;
    const QI_ARRAY(double) data;

    template <typename T> bool operator()(const T *const vin, T *rin) const {
        Eigen::Map<QI_ARRAY(T)>                 r(rin, data.rows());
        const Eigen::Map<const QI_ARRAYN(T, 3)> v(vin);

        const auto calc = model.spgr_signal(v);
        r               = data - calc;
        return true;
    }
};

struct HIFIMPRAGECost {
    const HIFIModel &model;
    const QI_ARRAY(double) data;

    template <typename T> bool operator()(const T *const vin, T *rin) const {
        Eigen::Map<QI_ARRAY(T)>                 r(rin, data.rows());
        const Eigen::Map<const QI_ARRAYN(T, 3)> v(vin);

        const auto calc = model.mprage_signal(v);
        r               = data - calc;
        return true;
    }
};

struct HIFIFit {
    static const bool Blocked = false;
    static const bool Indexed = false;
    using InputType           = double;
    using OutputType          = double;
    using RMSErrorType        = double;
    using FlagType            = int;
    using ModelType           = HIFIModel;
    HIFIModel model;

    int input_size(const int i) const {
        switch (i) {
        case 0:
            return model.spgr.size();
        case 1:
            return model.mprage.size();
        default:
            QI::Fail("Invalid input size = {}", i);
        }
    }
    int n_outputs() const { return 3; }

    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs,
                          HIFIModel::FixedArray const & /* Unused */,
                          HIFIModel::VaryingArray &    v,
                          HIFIModel::CovarArray *      cov,
                          RMSErrorType &               rmse,
                          std::vector<Eigen::ArrayXd> &residuals,
                          FlagType &                   iterations) const {
        double scale = std::max(inputs[0].maxCoeff(), inputs[1].maxCoeff());
        if (scale < std::numeric_limits<double>::epsilon()) {
            v << 0.0, 0.0, 0.0;
            rmse = 0.0;
            return {false, "Maximum data value was zero or less"};
        }
        const Eigen::ArrayXd spgr_data   = inputs[0] / scale;
        const Eigen::ArrayXd mprage_data = inputs[1] / scale;
        v << 10., 1., 1.; // PD, T1, B1
        ceres::Problem problem;
        using AutoSPGRType =
            ceres::AutoDiffCostFunction<HIFISPGRCost, ceres::DYNAMIC, HIFIModel::NV>;
        using AutoMPRAGEType =
            ceres::AutoDiffCostFunction<HIFIMPRAGECost, ceres::DYNAMIC, HIFIModel::NV>;
        auto *spgr_cost = new AutoSPGRType(new HIFISPGRCost{model, spgr_data}, model.spgr.size());
        auto *mprage_cost =
            new AutoMPRAGEType(new HIFIMPRAGECost{model, mprage_data}, model.mprage.size());
        problem.AddResidualBlock(spgr_cost, NULL, v.data());
        problem.AddResidualBlock(mprage_cost, NULL, v.data());
        for (int i = 0; i < 3; i++) {
            problem.SetParameterLowerBound(v.data(), i, model.bounds_lo[i]);
            problem.SetParameterUpperBound(v.data(), i, model.bounds_hi[i]);
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

        Eigen::ArrayXd const spgr_resid   = spgr_data - model.spgr_signal(v);
        Eigen::ArrayXd const mprage_resid = mprage_data - model.mprage_signal(v);
        if (residuals.size() > 0) {
            residuals[0] = spgr_resid * scale;
            residuals[1] = mprage_resid * scale;
        }
        double const var   = spgr_resid.square().sum() + mprage_resid.square().sum();
        int const    dsize = model.spgr.size() + model.mprage.size();
        if (cov) {
            QI::GetModelCovariance<ModelType>(problem, v, var / (dsize - ModelType::NV), cov);
        }
        rmse = sqrt(var / dsize);

        v[0] = v[0] * scale;
        return {true, ""};
    }
};

//******************************************************************************
// Main
//******************************************************************************
int despot1hifi_main(args::Subparser &parser) {
    args::Positional<std::string> spgr_path(parser, "SPGR_FILE", "Input SPGR file");
    args::Positional<std::string> mprage_path(parser, "MPRAGE_FILE", "Input MP-RAGE file");
    QI_COMMON_ARGS;
    args::ValueFlag<float> clamp(parser,
                                 "CLAMP",
                                 "Clamp output T1 values to this value",
                                 {'c', "clamp"},
                                 std::numeric_limits<float>::infinity());
    parser.Parse();

    QI::Log(verbose, "Reading sequence information");
    json input          = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    auto spgrSequence   = input.at("SPGR").get<QI::SPGRSequence>();
    auto mprageSequence = input.at("MPRAGE").get<QI::MPRAGESequence>();

    HIFIModel model{{}, spgrSequence, mprageSequence};
    if (simulate) {
        QI::SimulateModel<HIFIModel, true>(input,
                                           model,
                                           {},
                                           {QI::CheckPos(spgr_path), QI::CheckPos(mprage_path)},
                                           mask.Get(),
                                           verbose,
                                           simulate.Get(),
                                           subregion.Get());
    } else {
        HIFIFit hifi_fit{model};
        auto    fit_filter =
            QI::ModelFitFilter<HIFIFit>::New(&hifi_fit, verbose, covar, resids, subregion.Get());
        fit_filter->ReadInputs(
            {QI::CheckPos(spgr_path), QI::CheckPos(mprage_path)}, {}, mask.Get());
        fit_filter->Update();
        fit_filter->WriteOutputs(prefix.Get() + "HIFI_");
        QI::Log(verbose, "Finished.");
    }
    return EXIT_SUCCESS;
}
