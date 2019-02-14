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
#include "MPRAGESequence.h"
#include "Model.h"
#include "ModelFitFilter.h"
#include "OnePoolSignals.h"
#include "SPGRSequence.h"
#include "SequenceGroup.h"
#include "SimulateModel.h"
#include "Util.h"

using namespace std::literals;

struct HIFIModel {
    using SequenceType  = QI::SequenceGroup;
    using DataType      = double;
    using ParameterType = double;

    static constexpr int NV = 3;
    static constexpr int ND = 0;
    static constexpr int NF = 0;

    static std::array<const std::string, 3> varying_names;
    static std::array<const std::string, 0> fixed_names;
    static const QI_ARRAYN(double, 0) fixed_defaults;

    QI::SPGRSequence   spgr;
    QI::MPRAGESequence mprage;

    QI_ARRAYN(double, 3) bounds_lo = QI_ARRAYN(double, NV)::Constant(1.0e-4);
    QI_ARRAYN(double, 3)
    bounds_hi = QI_ARRAYN(double, NV)::Constant(std::numeric_limits<double>::infinity());

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
    auto spgr_signal(const Eigen::ArrayBase<Derived> &v, const QI_ARRAYN(double, NF) &
                     /* Unused */) const -> QI_ARRAY(typename Derived::Scalar) {
        return QI::SPGRSignal(v[0], v[1], v[2], spgr);
    }

    template <typename Derived>
    auto mprage_signal(const Eigen::ArrayBase<Derived> &v, const QI_ARRAYN(double, NF) &
                       /* Unused */) const -> QI_ARRAY(typename Derived::Scalar) {
        return QI::MPRAGESignal(v[0], v[1], v[2], mprage);
    }

    auto signals(const QI_ARRAYN(double, NV) & v, const QI_ARRAYN(double, NF) & f) const
        -> std::vector<QI_ARRAY(double)> {
        return {spgr_signal(v, f), mprage_signal(v, f)};
    }
};
std::array<const std::string, 3> HIFIModel::varying_names{{"PD"s, "T1"s, "B1"s}};
std::array<const std::string, 0> HIFIModel::fixed_names{{}};
const QI_ARRAYN(double, 0) HIFIModel::fixed_defaults{};

struct HIFISPGRCost {
    const HIFIModel &model;
    const QI_ARRAY(double) data;

    template <typename T> bool operator()(const T *const vin, T *rin) const {
        Eigen::Map<QI_ARRAY(T)>                 r(rin, data.rows());
        const Eigen::Map<const QI_ARRAYN(T, 3)> v(vin);
        QI_ARRAYN(double, 0) fixed;
        const auto calc = model.spgr_signal(v, fixed);
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
        QI_ARRAYN(double, 0) fixed;
        const auto calc = model.mprage_signal(v, fixed);
        r               = data - calc;
        return true;
    }
};

struct HIFIFit {
    static const bool Blocked = false;
    static const bool Indexed = false;
    using InputType           = double;
    using OutputType          = double;
    using ResidualType        = double;
    using FlagType            = int;
    using ModelType           = HIFIModel;
    HIFIModel model;

    int n_inputs() const { return 2; }
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
    int n_fixed() const { return 0; }
    int n_outputs() const { return 3; }

    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs,
                          const Eigen::ArrayXd & /* Unused */,
                          QI_ARRAYN(OutputType, HIFIModel::NV) & v,
                          ResidualType &               residual,
                          std::vector<Eigen::ArrayXd> &residuals,
                          FlagType &                   iterations) const {
        double scale = std::max(inputs[0].maxCoeff(), inputs[1].maxCoeff());
        if (scale < std::numeric_limits<double>::epsilon()) {
            v << 0.0, 0.0, 0.0;
            residual = 0.0;
            return std::make_tuple(false, "Maximum data value was zero or less");
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
            return std::make_tuple(false, summary.FullReport());
        }
        v[0]       = v[0] * scale;
        iterations = summary.iterations.size();
        residual   = summary.final_cost * scale;
        if (residuals.size() > 0) {
            std::vector<double> r_temp(model.spgr.size() + model.mprage.size());
            problem.Evaluate(ceres::Problem::EvaluateOptions(), NULL, &r_temp, NULL, NULL);
            for (int i = 0; i < model.spgr.size(); i++)
                residuals[0][i] = r_temp[i] * scale;
            for (int i = 0; i < model.mprage.size(); i++)
                residuals[1][i] = r_temp[i + model.spgr.size()] * scale;
        }
        return std::make_tuple(true, "");
    }
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser          parser("Calculates T1 and B1 maps from SPGR & IR-SPGR or MP-RAGE "
                                "data.\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> spgr_path(parser, "SPGR_FILE", "Input SPGR file");
    args::Positional<std::string> mprage_path(parser, "MPRAGE_FILE", "Input MP-RAGE file");
    QI_COMMON_ARGS;
    args::ValueFlag<float> clamp(parser,
                                 "CLAMP",
                                 "Clamp output T1 values to this value",
                                 {'c', "clamp"},
                                 std::numeric_limits<float>::infinity());
    QI::ParseArgs(parser, argc, argv, verbose, threads);

    QI::Log(verbose, "Reading sequence information");
    rapidjson::Document input = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    QI::SPGRSequence    spgrSequence(QI::GetMember(input, "SPGR"));
    QI::MPRAGESequence  mprageSequence(QI::GetMember(input, "MPRAGE"));
    HIFIModel           model{spgrSequence, mprageSequence};
    if (simulate) {
        QI::SimulateModel<HIFIModel, true>(input,
                                           model,
                                           {},
                                           {QI::CheckPos(spgr_path), QI::CheckPos(mprage_path)},
                                           verbose,
                                           simulate.Get());
    } else {
        HIFIFit hifi_fit{model};
        auto    fit_filter =
            QI::ModelFitFilter<HIFIFit>::New(&hifi_fit, verbose, resids, subregion.Get());
        fit_filter->ReadInputs(
            {QI::CheckPos(spgr_path), QI::CheckPos(mprage_path)}, {}, mask.Get());
        fit_filter->Update();
        fit_filter->WriteOutputs(prefix.Get() + "HIFI_");
        QI::Log(verbose, "Finished.");
    }
    return EXIT_SUCCESS;
}
