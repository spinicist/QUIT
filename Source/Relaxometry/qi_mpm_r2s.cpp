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
#include "Util.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

struct MPMModel {
    using SequenceType  = QI::MultiEchoSequence;
    using DataType      = double;
    using ParameterType = double;

    static constexpr int NV = 4;
    static constexpr int ND = 0;
    static constexpr int NF = 0;
    static constexpr int NI = 3;

    using VaryingArray = QI_ARRAYN(ParameterType, NV);
    using FixedArray   = QI_ARRAYN(ParameterType, NF);

    SequenceType &                    pdw_s, &t1w_s, &mtw_s;
    VaryingArray const                bounds_lo{1e-6, 1e-6, 1e-6, 1e-6};
    VaryingArray const                bounds_hi{1e4, 1e2, 1e2, 1e2}; // Signal values will be scaled
    std::array<std::string, NV> const varying_names{"R2s", "S0_PDw", "S0_T1w", "S0_MTw"};
    std::array<std::string, 0> const  fixed_names{};
    FixedArray const                  fixed_defaults{};

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

        r = data - model.pdw_signal(v);
        return true;
    }
};

struct T1wCost {
    MPMModel const &model;
    QI_ARRAY(double) const data;

    template <typename T> bool operator()(const T *const vin, T *rin) const {
        Eigen::Map<QI_ARRAYN(T, MPMModel::NV) const> const v(vin);

        Eigen::Map<QI_ARRAY(T)> r(rin, data.rows());

        r = data - model.t1w_signal(v);
        return true;
    }
};

struct MTwCost {
    MPMModel const &model;
    QI_ARRAY(double) const data;

    template <typename T> bool operator()(const T *const vin, T *rin) const {
        Eigen::Map<QI_ARRAYN(T, MPMModel::NV) const> const v(vin);

        Eigen::Map<QI_ARRAY(T)> r(rin, data.rows());

        r = data - model.mtw_signal(v);
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
                          RMSErrorType &               residual,
                          std::vector<Eigen::ArrayXd> &residuals,
                          FlagType &                   iterations) const {
        double scale = std::max({inputs[0].maxCoeff(), inputs[1].maxCoeff(), inputs[2].maxCoeff()});
        if (scale < std::numeric_limits<double>::epsilon()) {
            v        = ModelType::VaryingArray::Zero();
            residual = 0.0;
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
        residual   = summary.final_cost * scale;
        if (residuals.size() > 0) {
            residuals[0] = (pdw_data - model.pdw_signal(v)) * scale;
            residuals[1] = (t1w_data - model.t1w_signal(v)) * scale;
            residuals[2] = (mtw_data - model.mtw_signal(v)) * scale;
        }
        v.tail(3) = v.tail(3) * scale; // Multiply signals/proton densities back up
        return {true, ""};
    }
};

/*
 * Main
 */
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser(
        "Calculates R2* and S0 from PDw, T1w, MTw data.\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> pdw_path(parser, "PDw", "Input multi-echo PD-weighted file");
    args::Positional<std::string> t1w_path(parser, "T1w", "Input multi-echo T1-weighted file");
    args::Positional<std::string> mtw_path(parser, "MTw", "Input multi-echo MT-weighted file");

    args::HelpFlag       help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag           resids(parser, "RESIDS", "Write point residuals", {'r', "resids"});
    args::ValueFlag<int> threads(parser,
                                 "THREADS",
                                 "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"},
                                 QI::GetDefaultThreads());
    args::ValueFlag<std::string> prefix(
        parser, "OUTPREFIX", "Add a prefix to output filename", {'o', "out"});
    args::ValueFlag<std::string> mask_path(
        parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> subregion(
        parser,
        "SUBREGION",
        "Process voxels in a block from I,J,K with size SI,SJ,SK",
        {'s', "subregion"});
    args::ValueFlag<std::string> json_file(
        parser, "JSON", "Read JSON from file instead of stdin", {"json"});
    QI::ParseArgs(parser, argc, argv, verbose, threads);

    QI::CheckPos(pdw_path);
    QI::CheckPos(t1w_path);
    QI::CheckPos(mtw_path);

    QI::Log(verbose, "Reading sequence parameters");
    json                  doc = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);
    QI::MultiEchoSequence pdw_seq(doc["PDw"]), t1w_seq(doc["T1w"]), mtw_seq(doc["MTw"]);

    MPMModel model{pdw_seq, t1w_seq, mtw_seq};
    MPMFit   mpm_fit{model};
    auto fit_filter = QI::ModelFitFilter<MPMFit>::New(&mpm_fit, verbose, resids, subregion.Get());
    fit_filter->ReadInputs({pdw_path.Get(), t1w_path.Get(), mtw_path.Get()}, {}, mask_path.Get());
    fit_filter->Update();
    fit_filter->WriteOutputs(prefix.Get() + "MPM_");
    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
