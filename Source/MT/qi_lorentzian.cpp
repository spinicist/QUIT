/*
 *  qi_cestlorentz.cpp
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "ceres/ceres.h"
#include <Eigen/Core>
#include <iostream>

#include "Args.h"
#include "FitFunction.h"
#include "ImageIO.h"
#include "JSON.h"
#include "Model.h"
#include "ModelFitFilter.h"
#include "SequenceBase.h"
#include "SimulateModel.h"
#include "Util.h"

using namespace std::literals;

struct LorentzSequence : QI::SequenceBase {
    Eigen::ArrayXd sat_f0;
    QI_SEQUENCE_DECLARE(Lorentz);
    Eigen::Index size() const override { return sat_f0.rows(); }
};
LorentzSequence::LorentzSequence(const rapidjson::Value &json) {
    if (json.IsNull())
        QI_FAIL("Could not read sequence: " << name());
    sat_f0 = QI::ArrayFromJSON(json, "sat_f0");
}
rapidjson::Value LorentzSequence::toJSON(rapidjson::Document::AllocatorType &a) const {
    rapidjson::Value json(rapidjson::kObjectType);
    json.AddMember("sat_f0", QI::ArrayToJSON(sat_f0, a), a);
    return json;
}

using LorentzModel = QI::Model<4, 0, LorentzSequence>;
template <>
std::array<const std::string, 4> LorentzModel::varying_names{{"PD"s, "f0"s, "fwhm"s, "A"s}};
template <> std::array<const std::string, 0> LorentzModel::fixed_names{{}};
template <> const QI_ARRAYN(double, 0) LorentzModel::fixed_defaults{};
template <>
template <typename Derived>
auto LorentzModel::signal(const Eigen::ArrayBase<Derived> &v, const QI_ARRAYN(double, NF) &
                          /* Unused */) const -> QI_ARRAY(typename Derived::Scalar) {
    using T         = typename Derived::Scalar;
    const T &  PD   = v[0];
    const T &  f0   = v[1];
    const T &  fwhm = v[2];
    const T &  A    = v[3];
    const auto x    = (f0 - sequence.sat_f0) / (fwhm / 2.0);
    const auto L    = A / (1.0 + x.square());
    const auto s    = PD * (1.0 - L);
    return s;
}

struct LorentzFit : QI::FitFunction<LorentzModel> {
    using FitFunction::FitFunction;
    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs, const Eigen::ArrayXd &fixed,
                          QI_ARRAYN(OutputType, LorentzModel::NV) & p, ResidualType &      residual,
                          std::vector<Eigen::ArrayXd> & /* Unused */,
                          FlagType & /* Unused */) const override {
        const double scale = inputs[0].maxCoeff();
        if (scale < std::numeric_limits<double>::epsilon()) {
            p << 0.0, 0.0, 0.0, 0.0;
            residual = 0;
            return std::make_tuple(false, "Maximum data value was not positive");
        }
        const Eigen::ArrayXd &z_spec = inputs[0] / scale;

        using LCost     = QI::ModelCost<LorentzModel>;
        using AutoCost  = ceres::AutoDiffCostFunction<LCost, ceres::DYNAMIC, LorentzModel::NV>;
        auto *cost      = new LCost(model, fixed, z_spec);
        auto *auto_cost = new AutoCost(cost, this->model.sequence.size());
        //{"PD"s, "f0"s, "fwhm"s, "A"s}
        p << 2.0, 0.0, 0.5, 0.9;
        ceres::Problem problem;
        problem.AddResidualBlock(auto_cost, NULL, p.data());
        problem.SetParameterLowerBound(p.data(), 0, 1.e-6);
        problem.SetParameterUpperBound(p.data(), 0, 10.0);
        problem.SetParameterLowerBound(p.data(), 1, -10.0);
        problem.SetParameterUpperBound(p.data(), 1, 10.0);
        problem.SetParameterLowerBound(p.data(), 2, 1.e-6);
        problem.SetParameterUpperBound(p.data(), 2, 10.);
        problem.SetParameterLowerBound(p.data(), 3, 1.e-3);
        problem.SetParameterUpperBound(p.data(), 3, 10.0);
        ceres::Solver::Options options;
        options.function_tolerance  = 1e-5;
        options.gradient_tolerance  = 1e-6;
        options.parameter_tolerance = 1e-4;
        options.logging_type        = ceres::SILENT;
        ceres::Solver::Summary summary;
        options.max_num_iterations  = 50;
        options.function_tolerance  = 1e-5;
        options.gradient_tolerance  = 1e-6;
        options.parameter_tolerance = 1e-4;
        ceres::Solve(options, &problem, &summary);
        p[0] *= scale;
        residual = summary.final_cost;
        return std::make_tuple(true, "");
    }
};

int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Simple Lorentzian fitting.\nhttp://github.com/spinicist/QUIT");

    args::Positional<std::string> input_path(parser, "INPUT", "Input Z-spectrum file");
    args::HelpFlag                help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag           verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)",
                                 {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(parser, "PREFIX", "Add a prefix to output filenames",
                                        {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask",
                                      {'m', "mask"});
    args::ValueFlag<std::string> subregion(
        parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK",
        {'s', "subregion"});
    args::ValueFlag<std::string> seq_arg(parser, "FILE",
                                         "Read JSON input from file instead of stdin", {"file"});
    args::ValueFlag<float>       simulate(
        parser, "SIMULATE", "Simulate sequence instead of fitting model (argument is noise level)",
        {"simulate"}, 0.0);
    QI::ParseArgs(parser, argc, argv, verbose, threads);
    QI::CheckPos(input_path);
    QI_LOG(verbose, "Reading sequence information");
    rapidjson::Document input = seq_arg ? QI::ReadJSON(seq_arg.Get()) : QI::ReadJSON(std::cin);
    LorentzSequence     sequence(QI::GetMember(input, "Lorentz"));
    LorentzModel        model{sequence};
    if (simulate) {
        QI::SimulateModel<LorentzModel, false>(input, model, {}, {input_path.Get()}, verbose,
                                               simulate.Get());
    } else {
        LorentzFit fit{model};
        auto       fit_filter = itk::ModelFitFilter<LorentzFit>::New(&fit, verbose, false);
        fit_filter->SetInput(0, QI::ReadVectorImage(input_path.Get(), verbose));
        if (mask)
            fit_filter->SetMask(QI::ReadImage(mask.Get(), verbose));
        if (subregion)
            fit_filter->SetSubregion(QI::RegionArg(args::get(subregion)));
        fit_filter->Update();
        std::string outPrefix = outarg.Get() + "LTZ_";
        for (int i = 0; i < LorentzModel::NV; i++) {
            QI::WriteImage(fit_filter->GetOutput(i),
                           outPrefix + LorentzModel::varying_names.at(i) + QI::OutExt(), verbose);
        }
        QI::WriteImage(fit_filter->GetResidualOutput(), outPrefix + "residual" + QI::OutExt());
        QI_LOG(verbose, "Finished.");
    }
    return EXIT_SUCCESS;
}
