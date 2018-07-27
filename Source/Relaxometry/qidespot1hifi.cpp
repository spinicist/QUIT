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

#include <array>
#include <iostream>
#include <Eigen/Core>
#include "ceres/ceres.h"

#include "Model.h"
#include "FitFunction.h"
#include "ModelFitFilter.h"
#include "SimulateModel.h"
#include "SPGRSequence.h"
#include "MPRAGESequence.h"
#include "OnePoolSignals.h"
#include "Util.h"
#include "Args.h"
#include "ImageIO.h"

using namespace std::literals;

struct HIFIModel {
    using SequenceType  = QI::SequenceBase;
    using DataType      = double;
    using ParameterType = double;

    static const int NV = 3;
    static const int NF = 0;
    static std::array<const std::string, 3> varying_names;
    static std::array<const std::string, 0> fixed_names;
    static const QI_ARRAYN(double, 0) fixed_defaults;

    QI_ARRAYN(double, 3) bounds_lo = QI_ARRAYN(double, NV)::Constant(1.0e-4);
    QI_ARRAYN(double, 3) bounds_hi = QI_ARRAYN(double, NV)::Constant(std::numeric_limits<double>::infinity());

    template<typename Derived>
    auto spgr_signal(const Eigen::ArrayBase<Derived> &v,
                        const QI_ARRAYN(double, NF) &/* Unused */,
                        const QI::SPGRSequence *s) const -> QI_ARRAY(typename Derived::Scalar)
    {
        return QI::SPGRSignal(v[0], v[1], v[2], s);
    }

    template<typename Derived>
    auto mprage_signal(const Eigen::ArrayBase<Derived> &v,
                        const QI_ARRAYN(double, NF) &/* Unused */,
                        const QI::MPRAGESequence *s) const -> QI_ARRAY(typename Derived::Scalar)
    {
        return QI::MPRAGESignal(v[0], v[1], v[2], s);
    }

    template<typename Derived>
    auto signal(const Eigen::ArrayBase<Derived> &v,
                const QI_ARRAYN(double, NF) &f,
                const QI::SequenceBase *s) const -> QI_ARRAY(typename Derived::Scalar)
    {
        const QI::SPGRSequence *spgr = dynamic_cast<const QI::SPGRSequence *>(s);
        if (spgr) return spgr_signal(v, f, spgr);
        const QI::MPRAGESequence *mprage = dynamic_cast<const QI::MPRAGESequence *>(s);
        if (mprage) return mprage_signal(v, f, mprage);
        QI_FAIL("Given pointer was not to SPGRSequence or MPRAGESequence");
    }
};
std::array<const std::string, 3> HIFIModel::varying_names{{"PD"s, "T1"s, "B1"s}};
std::array<const std::string, 0> HIFIModel::fixed_names{{}};
const QI_ARRAYN(double, 0) HIFIModel::fixed_defaults{};

struct HIFISPGRCost {
    const HIFIModel &model;
    const QI::SPGRSequence *sequence;
    const QI_ARRAY(double) data;

    template<typename T>
    bool operator() (const T *const vin, T* rin) const {
        Eigen::Map<QI_ARRAY(T)> r(rin, data.rows());
        const Eigen::Map<const QI_ARRAYN(T, 3)> v(vin);
        QI_ARRAYN(double, 0) fixed;
        const auto calc = model.spgr_signal(v, fixed, sequence);
        r = data - calc;
        return true;
    }
};

struct HIFIMPRAGECost {
    const HIFIModel &model;
    const QI::MPRAGESequence *sequence;
    const QI_ARRAY(double) data;

    template<typename T>
    bool operator() (const T *const vin, T* rin) const {
        Eigen::Map<QI_ARRAY(T)> r(rin, data.rows());
        const Eigen::Map<const QI_ARRAYN(T, 3)> v(vin);
        QI_ARRAYN(double, 0) fixed;
        const auto calc = model.mprage_signal(v, fixed, sequence);
        r = data - calc;
        return true;
    }
};

struct HIFIFit {
    static const bool Blocked = false;
    static const bool Indexed = false;
    using InputType = double;
    using OutputType = double;
    using ResidualType = double;
    using FlagType = int;
    using ModelType = HIFIModel;
    QI::SPGRSequence *spgr_sequence;
    QI::MPRAGESequence *mprage_sequence;
    HIFIModel model;

    int n_inputs() const  { return 2; }
    int input_size(const int i) const {
        switch (i) {
        case 0: return spgr_sequence->size();
        case 1: return mprage_sequence->size();
        default: QI_FAIL("Invalid input " << i << " size requested");
        }
    }
    int n_fixed() const { return 0; }
    int n_outputs() const { return 3; }

    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs,
                          const Eigen::ArrayXd &/* Unused */, QI_ARRAYN(OutputType, HIFIModel::NV) &v,
                          ResidualType &residual, std::vector<Eigen::ArrayXd> &residuals, FlagType &iterations) const
    {
        double scale = std::max(inputs[0].maxCoeff(), inputs[1].maxCoeff());
        if (scale < std::numeric_limits<double>::epsilon()) {
            v << 0.0, 0.0, 0.0;
            residual = 0.0;
            return std::make_tuple(false, "Maximum data value was zero or less");
        }
        const Eigen::ArrayXd spgr_data = inputs[0] / scale;
        const Eigen::ArrayXd mprage_data = inputs[1] / scale;
        v << 10., 1., 1.; // PD, T1, B1
        ceres::Problem problem;
        using AutoSPGRType = ceres::AutoDiffCostFunction<HIFISPGRCost, ceres::DYNAMIC, HIFIModel::NV>;
        using AutoMPRAGEType = ceres::AutoDiffCostFunction<HIFIMPRAGECost, ceres::DYNAMIC, HIFIModel::NV>;
        auto *spgr_cost = new AutoSPGRType(new HIFISPGRCost{model, spgr_sequence, spgr_data}, spgr_sequence->size());
        auto *mprage_cost = new AutoMPRAGEType(new HIFIMPRAGECost{model, mprage_sequence, mprage_data}, mprage_sequence->size());
        problem.AddResidualBlock(spgr_cost, NULL, v.data());
        problem.AddResidualBlock(mprage_cost, NULL, v.data());
        for (int i = 0; i < 3; i++) {
            problem.SetParameterLowerBound(v.data(), i, model.bounds_lo[i]);
            problem.SetParameterUpperBound(v.data(), i, model.bounds_hi[i]);
        }
        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
        options.max_num_iterations = 50;
        options.function_tolerance = 1e-5;
        options.gradient_tolerance = 1e-6;
        options.parameter_tolerance = 1e-4;
        options.logging_type = ceres::SILENT;
        ceres::Solve(options, &problem, &summary);
        if (!summary.IsSolutionUsable()) {
            return std::make_tuple(false, summary.FullReport());
        }
        v[0] = v[0] * scale;
        iterations = summary.iterations.size();
        residual = summary.final_cost * scale;
        if (residuals.size() > 0) {
            std::vector<double> r_temp(spgr_sequence->size() + mprage_sequence->size());
            problem.Evaluate(ceres::Problem::EvaluateOptions(), NULL, &r_temp, NULL, NULL);
            for (int i = 0; i < spgr_sequence->size(); i++)
                residuals[0][i] = r_temp[i] * scale;
            for (int i = 0; i < mprage_sequence->size(); i++)
                residuals[1][i] = r_temp[i + spgr_sequence->size()] * scale;
        }
        return std::make_tuple(true, "");
    }
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates T1 and B1 maps from SPGR & IR-SPGR or MP-RAGE data.\nhttp://github.com/spinicist/QUIT");
    
    args::Positional<std::string> spgr_path(parser, "SPGR_FILE", "Input SPGR file");
    args::Positional<std::string> mprage_path(parser, "MPRAGE_FILE", "Input MP-RAGE file");
    
    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag     all_resids(parser, "ALL RESIDUALS", "Output individual residuals in addition to the Sum-of-Squares", {'r',"resids"});
    args::ValueFlag<float> clamp(parser, "CLAMP", "Clamp output T1 values to this value", {'c', "clamp"}, std::numeric_limits<float>::infinity());
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> subregion(parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    args::Flag resids(parser, "RESIDS", "Write out residuals for each data-point", {'r', "resids"});
    args::ValueFlag<std::string> seq_arg(parser, "FILE", "Read JSON input from file instead of stdin", {"file"});
    args::ValueFlag<float> simulate(parser, "SIMULATE", "Simulate sequence instead of fitting model (argument is noise level)", {"simulate"}, 0.0);
    QI::ParseArgs(parser, argc, argv, verbose, threads);

    QI_LOG(verbose, "Reading sequence information");
    rapidjson::Document input = seq_arg ? QI::ReadJSON(seq_arg.Get()) : QI::ReadJSON(std::cin);
    QI::SPGRSequence spgrSequence(QI::GetMember(input, "SPGR"));
    QI::MPRAGESequence mprageSequence(QI::GetMember(input, "MPRAGE"));

    if (simulate) {
        HIFIModel model;
        QI::SimulateModel<HIFIModel, false>(input, model, {&spgrSequence, &mprageSequence}, {}, {QI::CheckPos(spgr_path), QI::CheckPos(mprage_path)}, verbose, simulate.Get());
    } else {
        HIFIFit hifi_fit;
        hifi_fit.spgr_sequence = &spgrSequence;
        hifi_fit.mprage_sequence = &mprageSequence;
        auto fit_filter = itk::ModelFitFilter<HIFIFit>::New();
        fit_filter->SetVerbose(verbose);
        fit_filter->SetFitFunction(&hifi_fit);
        fit_filter->SetOutputAllResiduals(resids);
        fit_filter->SetInput(0, QI::ReadVectorImage(QI::CheckPos(spgr_path), verbose));
        fit_filter->SetInput(1, QI::ReadVectorImage(QI::CheckPos(mprage_path), verbose));
        if (mask) fit_filter->SetMask(QI::ReadImage(mask.Get(), verbose));
        if (subregion) fit_filter->SetSubregion(QI::RegionArg(args::get(subregion)));
        QI_LOG(verbose, "Processing");
        if (verbose) {
            auto monitor = QI::GenericMonitor::New();
            fit_filter->AddObserver(itk::ProgressEvent(), monitor);
        }
        fit_filter->Update();
        QI_LOG(verbose, "Elapsed time was " << fit_filter->GetTotalTime() << "s\n" <<
                        "Writing results files.");
        std::string outPrefix = outarg.Get() + "HIFI_";
        for (int i = 0; i < hifi_fit.model.NV; i++) {
            QI::WriteImage(fit_filter->GetOutput(i), outPrefix + hifi_fit.model.varying_names.at(i) + QI::OutExt());
        }
        QI::WriteImage(fit_filter->GetResidualOutput(), outPrefix + "residual" + QI::OutExt());
        if (resids) {
            QI::WriteVectorImage(fit_filter->GetResidualsOutput(0), outPrefix + "all_residuals" + QI::OutExt());
        }
        QI_LOG(verbose, "Finished." );
    }
    return EXIT_SUCCESS;
}
