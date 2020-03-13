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

struct JSRModel : QI::Model<double, double, 4, 1, 2> {
    // Sequence paramter structs
    QI::SPGREchoSequence &  spgr;
    QI::SSFPFiniteSequence &ssfp;

    // Fitting start point and bounds
    // The PD will be scaled by the fitting function to keep parameters roughly the same magnitude
    // The off-resonance fit is done on the accrued angle and then converted to frequency after
    VaryingArray const start{7., 1., 0.05, 0.};
    VaryingArray const bounds_lo{1, 1e-3, 1e-3, -2 * M_PI};
    VaryingArray const bounds_hi{15, 5, 3, 2 * M_PI};

    std::array<std::string, NV> const varying_names{"PD", "T1", "T2", "df0"};
    std::array<std::string, NF> const fixed_names{"B1"};
    // If fixed parameters not supplied, use these default values
    FixedArray const fixed_defaults{1.0};

    // Signal functions. These have to be templated to allow automatic differentiation within Ceres
    template <typename Derived>
    auto spgr_signal(Eigen::ArrayBase<Derived> const &v, FixedArray const &f) const
        -> QI_ARRAY(typename Derived::Scalar) {
        // Get the underlying datatype of the passed-in Eigen Array (usually double or a Ceres Jet)
        using T = typename Derived::Scalar;
        // Now get the individual parameters by name for clarity
        T const &PD = v[0];
        T const &T1 = v[1];
        T const &T2 = v[2];

        // Fixed parameters will never be Ceres Jets, so can use a raw type
        double const &B1             = f[0];
        QI_ARRAY(double) const alpha = spgr.FA * B1;

        // Anything expression that involves a varying parameter must be of type T
        // Conversely, don't try to initialise a double/T from a T/double (fun error messages)
        T const E1 = exp(-spgr.TR / T1);
        T const Ee = exp(-spgr.TE / T2);

        // The final signal will be an array of T, but the length is dependant on the sequence
        // parameters, so it is not a VaryingArray which has fixed length
        QI_ARRAY(T) const signal = PD * Ee * sin(alpha) * (1.0 - E1) / (1.0 - E1 * cos(alpha));
        return signal;
    }

    template <typename Derived>
    auto ssfp_signal(Eigen::ArrayBase<Derived> const &v, FixedArray const &f) const
        -> QI_ARRAY(typename Derived::Scalar) {
        using T      = typename Derived::Scalar;
        T const &PD  = v[0];
        T const &T1  = v[1];
        T const &T2  = v[2];
        T const &psi = v[3];

        double const &B1             = f[0];
        QI_ARRAY(double) const alpha = ssfp.FA * B1;

        // Crooijman / Bieri correction
        T const T_rfe = (0.68 - 0.125 * (1.0 + ssfp.Trf / ssfp.TR) * T2 / T1) * ssfp.Trf;
        T const TRc   = ssfp.TR - T_rfe;

        // This is a useful trick to get debugging information only when Ceres is evaluating the
        // cost function at the end of an iteration. Printing Jet objects on every evaluation gives
        // too much output
        // if constexpr (std::is_floating_point<T>::value) {
        //     QI_DB(B1)
        //     QI_DBVEC(ssfp.FA)
        //     QI_DBVEC(alpha)
        // }

        T const E1 = exp(-ssfp.TR / T1);
        T const E2 = exp(-TRc / T2);
        T const Ee = exp(-TRc / (2.0 * T2));

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
        QI_ARRAY(T) const s = sqrt(re_m.square() + im_m.square());
        return s;
    }

    auto signals(VaryingArray const &v, FixedArray const &f) const
        -> std::vector<QI_ARRAY(double)> {
        return {spgr_signal(v, f), ssfp_signal(v, f)};
    }
};

// Cost functors. These need to calculate the residuals
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

// Fit function structure. This is what actually runs the fitting/optimisation
struct JSRFit {
    // Boilerplate information required by ModelFitFilter
    static const bool Blocked = false; // = input is in blocks and outputs have multiple entries
    static const bool Indexed = false; // = the voxel index will be passed to the fit
    using RMSErrorType        = double;
    using FlagType            = int; // Almost always the number of iterations

    using ModelType = JSRModel;
    ModelType model;
    int       n_psi;

    // Have to tell the ModelFitFilter how many volumes we expect in each input
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

    // This has to match the function signature that will be called in ModelFitFilter (which depends
    // on Blocked/Indexed. The return type is a simple struct indicating success, and on failure
    // also the reason for failure
    QI::FitReturnType
    fit(std::vector<Eigen::ArrayXd> const &inputs,       // Input: signal data
        ModelType::FixedArray const &      fixed,        // Input: Fixed parameters
        ModelType::VaryingArray &          best_varying, // Output: Varying parameters
        ModelType::CovarArray *            covar,
        RMSErrorType &                     rmse,      // Output: root-mean-square error
        std::vector<Eigen::ArrayXd> &      residuals, // Optional output: point residuals
        FlagType &                         iterations /* Usually iterations */) const {
        // First scale down the raw data so that PD will be roughly the same magnitude as other
        // parameters This is important for numerical stability in the optimiser

        double scale = std::max({inputs[0].maxCoeff(), inputs[1].maxCoeff()});
        if (scale < std::numeric_limits<double>::epsilon()) {
            best_varying = ModelType::VaryingArray::Zero();
            rmse         = 0.0;
            return {false, "Maximum data value was zero or less"};
        }
        Eigen::ArrayXd const spgr_data = inputs[0] / scale;
        Eigen::ArrayXd const ssfp_data = inputs[1] / scale;

        // Setup Ceres
        ceres::Problem problem;
        using AutoSPGRType = ceres::AutoDiffCostFunction<SPGRCost, ceres::DYNAMIC, ModelType::NV>;
        using AutoSSFPType = ceres::AutoDiffCostFunction<SSFPCost, ceres::DYNAMIC, ModelType::NV>;
        auto *spgr_cost =
            new AutoSPGRType(new SPGRCost{model, fixed, spgr_data}, model.spgr.size());
        auto *ssfp_cost =
            new AutoSSFPType(new SSFPCost{model, fixed, ssfp_data}, model.ssfp.size());
        ceres::LossFunction *loss = new ceres::HuberLoss(1.0); // Don't know if this helps
        // This is where the parameters and cost functions actually get added to Ceres
        ModelType::VaryingArray varying;
        problem.AddResidualBlock(spgr_cost, loss, varying.data());
        problem.AddResidualBlock(ssfp_cost, loss, varying.data());

        // Set up parameter bounds
        for (int i = 0; i < ModelType::NV; i++) {
            problem.SetParameterLowerBound(varying.data(), i, model.bounds_lo[i]);
            problem.SetParameterUpperBound(varying.data(), i, model.bounds_hi[i]);
        }

        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
        options.max_num_iterations  = 50;
        options.function_tolerance  = 1e-6;
        options.gradient_tolerance  = 1e-7;
        options.parameter_tolerance = 1e-5;
        options.logging_type        = ceres::SILENT;

        // We need to do 2 starts for JSR in case off-resonance is very high
        double       best_cost = std::numeric_limits<double>::max();
        double const psi_step  = (n_psi % 2) ? 2 * M_PI / (n_psi - 1) : 2 * M_PI / (n_psi);
        double       psi       = (n_psi == 1) ? 0 : -M_PI;
        for (int p = 0; p < n_psi; p++, psi += psi_step) {
            varying    = model.start;
            varying[3] = psi;
            ceres::Solve(options, &problem, &summary);
            if (!summary.IsSolutionUsable()) {
                return {false, summary.FullReport()};
            }
            if (summary.final_cost < best_cost) {
                iterations   = summary.iterations.size();
                best_varying = varying;
                best_cost    = summary.final_cost;
            }
        }
        Eigen::ArrayXd const spgr_residual = (spgr_data - model.spgr_signal(best_varying, fixed));
        Eigen::ArrayXd const ssfp_residual = (ssfp_data - model.ssfp_signal(best_varying, fixed));
        if (residuals.size() > 0) {
            residuals[0] = spgr_residual * scale;
            residuals[1] = ssfp_residual * scale;
        }
        double const var   = spgr_residual.square().sum() + ssfp_residual.square().sum();
        int const    dsize = model.spgr.size() + model.ssfp.size();
        rmse               = sqrt(spgr_residual.square().mean() + ssfp_residual.square().mean());
        if (covar) {
            varying = best_varying;
            QI::GetModelCovariance<JSRModel>(problem, varying, var / (dsize - JSRModel::NV), covar);
        }
        best_varying[0] *= scale; // Multiply signals/proton density back up
        // Wrap and convert to frequency
        best_varying[3] =
            (std::fmod(best_varying[3] + 3 * M_PI, 2 * M_PI) - M_PI) / (2 * M_PI * model.ssfp.TR);
        return {true, ""};
    }
};

/*
 * Main
 */
int jsr_main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser          parser("Calculates T1/T2 from simultaneous fit to SPGR/SSFP "
                                "data.\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> spgr_path(parser, "SPGR", "Input SPGR file");
    args::Positional<std::string> ssfp_path(parser, "SSFP", "Input SSFP file");

    QI_COMMON_ARGS;

    args::ValueFlag<std::string> b1_path(parser, "B1", "Path to B1 map", {'b', "B1"});
    args::ValueFlag<int>         npsi(
        parser, "N PSI", "Number of starts for psi/off-resonance, default 2", {'p', "npsi"}, 2);

    QI::ParseArgs(parser, argc, argv, verbose, threads);

    QI::CheckPos(spgr_path);
    QI::CheckPos(ssfp_path);

    QI::Log(verbose, "Reading sequence parameters");
    json doc = json_file ? QI::ReadJSON(json_file.Get()) : QI::ReadJSON(std::cin);

    QI::SPGREchoSequence   spgr_seq(doc["SPGR"]);
    QI::SSFPFiniteSequence ssfp_seq(doc["SSFP"]);

    JSRModel model{{}, spgr_seq, ssfp_seq};
    JSRFit   jsr_fit{model, npsi.Get()};
    auto     fit_filter =
        QI::ModelFitFilter<JSRFit>::New(&jsr_fit, verbose, covar, resids, subregion.Get());
    fit_filter->ReadInputs({spgr_path.Get(), ssfp_path.Get()}, {b1_path.Get()}, mask.Get());
    fit_filter->Update();
    fit_filter->WriteOutputs(prefix.Get() + "JSR_");
    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
