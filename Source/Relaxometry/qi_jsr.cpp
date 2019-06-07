/*
 *  qi_jsr.cpp
 *
 *  Copyright (c) 2019 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <Eigen/Core>
#include <type_traits>

// #define QI_DEBUG_BUILD 1

#include "Args.h"
#include "ImageIO.h"
#include "Macro.h"
#include "Model.h"
#include "ModelFitFilter.h"
#include "SPGRSequence.h"
#include "SSFPSequence.h"
#include "Util.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

struct JSRModel {
    using DataType      = double;
    using ParameterType = double;

    static constexpr int NV = 4;
    static constexpr int ND = 0;
    static constexpr int NF = 1;
    static constexpr int NI = 2;

    using VaryingArray = QI_ARRAYN(ParameterType, NV);
    using FixedArray   = QI_ARRAYN(ParameterType, NF);

    QI::SPGREchoSequence &spgr;
    QI::SSFPSequence &    ssfp;

    VaryingArray const                start{10., 1., 0.1, 0.};
    VaryingArray const                bounds_lo{1e-6, 1e-3, 1e-3, -1.0 / ssfp.TR};
    VaryingArray const                bounds_hi{1e2, 5, 3, 1.0 / ssfp.TR};
    std::array<std::string, NV> const varying_names{"PD", "T1", "T2", "df0"};
    std::array<std::string, NF> const fixed_names{"B1"};
    FixedArray const                  fixed_defaults{1.0};

    template <typename Derived>
    auto spgr_signal(Eigen::ArrayBase<Derived> const &v, FixedArray const &f) const
        -> QI_ARRAY(typename Derived::Scalar) {
        using T     = typename Derived::Scalar;
        T const &PD = v[0];
        T const &T1 = v[1];
        T const &T2 = v[2];

        double const &B1             = f[0];
        QI_ARRAY(double) const alpha = spgr.FA * B1;

        T const E1 = exp(-spgr.TR / T1);
        T const Ee = exp(-spgr.TE / T2);

        QI_ARRAY(T) const signal = PD * Ee * sin(alpha) * (1.0 - E1) / (1.0 - E1 * cos(alpha));
        return signal;
    }

    template <typename Derived>
    auto ssfp_signal(Eigen::ArrayBase<Derived> const &v, FixedArray const &f) const
        -> QI_ARRAY(typename Derived::Scalar) {
        using T     = typename Derived::Scalar;
        T const &PD = v[0];
        T const &T1 = v[1];
        T const &T2 = v[2];
        T const &f0 = v[3];

        double const &B1             = f[0];
        QI_ARRAY(double) const alpha = ssfp.FA * B1;

        T const E1  = exp(-ssfp.TR / T1);
        T const E2  = exp(-ssfp.TR / T2);
        T const Ee  = exp(-ssfp.TR / (2.0 * T2));
        T const psi = 2. * M_PI * f0 * ssfp.TR;

        QI_ARRAY(T) const d = (1. - E1 * E2 * E2 - (E1 - E2 * E2) * cos(alpha));
        QI_ARRAY(T) const G = -PD * Ee * (1. - E1) * sin(alpha) / d;
        QI_ARRAY(T) const b = E2 * (1. - E1) * (1. + cos(alpha)) / d;

        QI_ARRAY(T) const theta  = ssfp.PhaseInc + psi;
        QI_ARRAY(T) const cos_th = cos(theta);
        QI_ARRAY(T) const sin_th = sin(theta);
        T const cos_psi          = cos(psi);
        T const sin_psi          = sin(psi);

        QI_ARRAY(T)
        const re_m =
            (cos_psi - E2 * (cos_th * cos_psi - sin_th * sin_psi)) * G / (1.0 - b * cos_th);
        QI_ARRAY(T)
        const im_m =
            (sin_psi - E2 * (cos_th * sin_psi + sin_th * cos_psi)) * G / (1.0 - b * cos_th);
        QI_ARRAY(T) const signal = sqrt(re_m.square() + im_m.square());
        return signal;
    }

    auto signals(VaryingArray const &v, FixedArray const &f) const
        -> std::vector<QI_ARRAY(double)> {
        return {spgr_signal(v, f), ssfp_signal(v, f)};
    }
};

struct SPGRCost {
    JSRModel const &     model;
    JSRModel::FixedArray fixed;
    QI_ARRAY(double) const data;

    template <typename T> bool operator()(T const *const vin, T *rin) const {
        Eigen::Map<QI_ARRAYN(T, JSRModel::NV) const> const varying(vin);

        Eigen::Map<QI_ARRAY(T)> residuals(rin, data.rows());
        residuals = data - model.spgr_signal(varying, fixed);
        return true;
    }
};

struct SSFPCost {
    JSRModel const &     model;
    JSRModel::FixedArray fixed;
    QI_ARRAY(double) const data;

    template <typename T> bool operator()(T const *const vin, T *rin) const {
        Eigen::Map<QI_ARRAYN(T, JSRModel::NV) const> const varying(vin);

        Eigen::Map<QI_ARRAY(T)> residuals(rin, data.rows());
        residuals = data - model.ssfp_signal(varying, fixed);
        return true;
    }
};

struct JSRFit {
    static const bool Blocked = false;
    static const bool Indexed = false;
    using InputType           = double;
    using OutputType          = double;
    using ResidualType        = double;
    using FlagType            = int;
    using ModelType           = JSRModel;
    ModelType model;

    int input_size(const int i) const {
        switch (i) {
        case 0:
            return model.spgr.size();
        case 1:
            return model.ssfp.size();
        default:
            QI::Fail("Invalid input size = {}", i);
        }
    }
    int n_outputs() const { return model.NV; }

    QI::FitReturnType fit(std::vector<Eigen::ArrayXd> const &inputs,
                          Eigen::ArrayXd const &             fixed,
                          ModelType::VaryingArray &          best_varying,
                          ResidualType &                     residual,
                          std::vector<Eigen::ArrayXd> &      residuals,
                          FlagType &                         iterations) const {
        double scale = std::max({inputs[0].maxCoeff(), inputs[1].maxCoeff()});
        if (scale < std::numeric_limits<double>::epsilon()) {
            best_varying = ModelType::VaryingArray::Zero();
            residual     = 0.0;
            return {false, "Maximum data value was zero or less"};
        }
        Eigen::ArrayXd const spgr_data = inputs[0] / scale;
        Eigen::ArrayXd const ssfp_data = inputs[1] / scale;

        ceres::Problem problem;
        using AutoSPGRType = ceres::AutoDiffCostFunction<SPGRCost, ceres::DYNAMIC, ModelType::NV>;
        using AutoSSFPType = ceres::AutoDiffCostFunction<SSFPCost, ceres::DYNAMIC, ModelType::NV>;
        auto *spgr_cost =
            new AutoSPGRType(new SPGRCost{model, fixed, spgr_data}, model.spgr.size());
        auto *ssfp_cost =
            new AutoSSFPType(new SSFPCost{model, fixed, ssfp_data}, model.ssfp.size());
        ceres::LossFunction *   loss = new ceres::HuberLoss(1.0);
        ModelType::VaryingArray varying;
        problem.AddResidualBlock(spgr_cost, loss, varying.data());
        problem.AddResidualBlock(ssfp_cost, loss, varying.data());
        for (int i = 0; i < ModelType::NV; i++) {
            problem.SetParameterLowerBound(varying.data(), i, model.bounds_lo[i]);
            problem.SetParameterUpperBound(varying.data(), i, model.bounds_hi[i]);
        }
        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
        options.max_num_iterations  = 50;
        options.function_tolerance  = 1e-5;
        options.gradient_tolerance  = 1e-6;
        options.parameter_tolerance = 1e-4;
        options.logging_type        = ceres::SILENT;

        double best_cost = std::numeric_limits<double>::max();
        for (double f0_start : {0., -0.5 / model.ssfp.TR, 0.5 / model.ssfp.TR}) {
            varying    = model.start;
            varying[3] = f0_start;
            ceres::Solve(options, &problem, &summary);
            if (!summary.IsSolutionUsable()) {
                return {false, summary.FullReport()};
            }
            if (summary.final_cost < best_cost) {
                iterations   = summary.iterations.size();
                residual     = summary.final_cost * scale;
                best_varying = varying;
                best_cost    = summary.final_cost;
            }
        }
        if (residuals.size() > 0) {
            residuals[0] = (spgr_data - model.spgr_signal(best_varying, fixed)) * scale;
            residuals[1] = (ssfp_data - model.ssfp_signal(best_varying, fixed)) * scale;
        }
        best_varying[0] *= scale; // Multiply signals/proton density back up
        return {true, ""};
    }
};

/*
 * Main
 */
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser          parser("Calculates T1/T2 from simultaneous fit to SPGR/SSFP "
                                "data.\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> spgr_path(parser, "SPGR", "Input SPGR file");
    args::Positional<std::string> ssfp_path(parser, "SSFP", "Input SSFP file");

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
    args::ValueFlag<std::string> b1_path(parser, "B1", "Path to B1 map", {'b', "B1"});
    QI::ParseArgs(parser, argc, argv, verbose, threads);

    QI::CheckPos(spgr_path);
    QI::CheckPos(ssfp_path);

    QI::Log(verbose, "Reading sequence parameters");
    json doc = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);

    QI::SPGREchoSequence spgr_seq(doc["SPGR"]);
    QI::SSFPSequence     ssfp_seq(doc["SSFP"]);

    JSRModel model{spgr_seq, ssfp_seq};
    JSRFit   jsr_fit{model};
    auto fit_filter = QI::ModelFitFilter<JSRFit>::New(&jsr_fit, verbose, resids, subregion.Get());
    fit_filter->ReadInputs({spgr_path.Get(), ssfp_path.Get()}, {}, mask_path.Get());
    fit_filter->Update();
    fit_filter->WriteOutputs(prefix.Get() + "JSR_");
    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
