/*
 *  qmcdespot.cpp
 *
 *  Created by Tobias Wood on 03/06/2015.
 *  Copyright (c) 2015 Tobias Wood.
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
#include "SSFPSequence.h"
#include "TwoPoolModel.h"
#include "ThreePoolModel.h"
#include "Util.h"
#include "Args.h"
#include "ImageIO.h"
#include "RegionContraction.h"

template<typename Model>
struct MCDSRCFunctor {
    const Eigen::ArrayXd data, weights;
    const QI_ARRAYN(double, Model::NF) fixed;
    const Model &model;

    MCDSRCFunctor(const Model &m, const QI_ARRAYN(double, Model::NF) &f,
                  const Eigen::ArrayXd &d, const Eigen::ArrayXd &w) :
        data(d), weights(w), fixed(f), model(m)
    {
        assert(static_cast<size_t>(data.rows()) == model.sequence.size());
    }

    int inputs() const { return Model::NV; }
    int values() const { return model.sequence.size(); }

    bool constraint(const QI_ARRAYN(double, Model::NV) &varying) const {
        return model.valid(varying);
    }

    Eigen::ArrayXd residuals(const QI_ARRAYN(double, Model::NV) &varying) const {
        return data - model.signal(varying, fixed);
    }

    double operator()(const QI_ARRAYN(double, Model::NV) &varying) const {
        return (residuals(varying) * weights).square().sum();
    }
};

template<typename Model>
struct SRCFit  {
    static const bool Blocked = false;
    static const bool Indexed = false;
    using InputType = double;
    using OutputType = double;
    using ResidualType = double;
    using FlagType = int;
    using ModelType = Model;
    Model &model;

    int n_inputs() const { return model.sequence.count(); }
    int input_size(const int i) const { return model.sequence.at(i)->size(); }
    int n_fixed() const { return Model::NF; }
    int n_outputs() const { return Model::NV; }
    
    int max_iterations = 5;
    size_t src_samples = 5000, src_retain = 50;
    bool src_gauss = true;

    QI::FitReturnType fit(const std::vector<Eigen::ArrayXd> &inputs,
                          const Eigen::ArrayXd &fixed, QI_ARRAYN(OutputType, Model::NV) &v,
                          ResidualType &residual, std::vector<Eigen::ArrayXd> &residuals, FlagType &iterations) const
    {
        Eigen::ArrayXd data(model.sequence.size());
        int dataIndex = 0;
        for (size_t i = 0; i < inputs.size(); i++) {
            if (model.scale_to_mean) {
                data.segment(dataIndex, inputs[i].rows()) = inputs[i] / inputs[i].mean();
            } else {
                data.segment(dataIndex, inputs[i].rows()) = inputs[i];
            }
            dataIndex += inputs[i].rows();
        }
        QI_ARRAYN(double, Model::NV) thresh = QI_ARRAYN(double, Model::NV)::Constant(0.05);
        const double &f0 = fixed[0];
        Eigen::ArrayXd weights = model.sequence.weights(f0);
        using Functor = MCDSRCFunctor<Model>;
        Functor func(model, fixed, data, weights);
        QI::RegionContraction<Functor> rc(func, model.bounds_lo, model.bounds_hi, thresh, src_samples, src_retain, max_iterations, 0.02, src_gauss, false);
        if (!rc.optimise(v)) {
            return std::make_tuple(false, "Region contraction failed");
        }
        auto r = func.residuals(v);
        residual = sqrt(r.square().sum() / r.rows());
        if (residuals.size() > 0) {
            int index = 0;
            for (size_t i = 0; i < model.sequence.count(); i++) {
                residuals[i] = r.segment(index, model.sequence.at(i)->size());
                index += model.sequence.at(i)->size();
            }
        }
        iterations = rc.contractions();
        return std::make_tuple(true, "");
    }
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates MWF & other parameter maps from mcDESPOT data\n"
                                "All times (e.g. T1, TR) are in SECONDS. All angles are in degrees.\n"
                                "http://github.com/spinicist/QUIT");
    args::PositionalList<std::string> input_paths(parser, "INPUT FILES", "Input image files");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> f0(parser, "f0", "f0 map (Hertz)", {'f', "f0"});
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio)", {'b', "B1"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> subregion(parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    args::Flag resids(parser, "RESIDS", "Write out residuals for each data-point", {'r', "resids"});
    args::ValueFlag<int> modelarg(parser, "MODEL", "Select model to fit - 2/3, default 3", {'M', "model"}, 3);
    args::Flag scale(parser, "SCALE", "Normalize signals to mean (a good idea)", {'S', "scale"});
    args::Flag use_src(parser, "SRC", "Use flat prior (stochastic region contraction), not gaussian", {"SRC"});
    args::ValueFlag<int> its(parser, "ITERS", "Max iterations, default 4", {'i',"its"}, 4);
    args::Flag bounds(parser, "BOUNDS", "Specify bounds in input", {"bounds"});
    args::ValueFlag<std::string> seq_arg(parser, "FILE", "Read JSON input from file instead of stdin", {"file"});
    args::ValueFlag<float> simulate(parser, "SIMULATE", "Simulate sequence instead of fitting model (argument is noise level)", {"simulate"}, 0.0);
    QI::ParseArgs(parser, argc, argv, verbose, threads);
    QI::CheckList(input_paths);

    QI_LOG(verbose, "Reading sequences");
    rapidjson::Document input = seq_arg ? QI::ReadJSON(seq_arg.Get()) : QI::ReadJSON(std::cin);
    QI::SequenceGroup sequences(QI::GetMember(input, "Sequences"));

    auto process = [&](auto model, const std::string &prefix) {
        if (simulate) {
            QI::SimulateModel<decltype(model), true>(input, model, {f0.Get(), B1.Get()}, input_paths.Get(), verbose, simulate.Get());
        } else {
            using FitType = SRCFit<decltype(model)>;
            FitType src{model};
            src.src_gauss = !use_src;
            if (bounds) {
                src.model.bounds_lo = QI::ArrayFromJSON(input, "lower_bounds");
                src.model.bounds_hi = QI::ArrayFromJSON(input, "upper_bounds");
            }
            QI_LOG(verbose, "Low bounds: " << src.model.bounds_lo.transpose() <<
                            "\nHigh bounds: " << src.model.bounds_hi.transpose());

            auto fit_filter = itk::ModelFitFilter<FitType>::New(&src);
            fit_filter->SetVerbose(verbose);
            fit_filter->SetOutputAllResiduals(resids);
            for (size_t i = 0; i < input_paths.Get().size(); i++) {
                fit_filter->SetInput(i, QI::ReadVectorImage(input_paths.Get().at(i), verbose));
            }
            if (f0) fit_filter->SetFixed(0, QI::ReadImage(f0.Get(), verbose));
            if (B1) fit_filter->SetFixed(1, QI::ReadImage(B1.Get(), verbose));
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
            std::string outPrefix = outarg.Get() + prefix;
            for (int i = 0; i < src.n_outputs(); i++) {
                QI::WriteImage(fit_filter->GetOutput(i), outPrefix + src.model.varying_names.at(i) + QI::OutExt());
            }
            QI::WriteImage(fit_filter->GetResidualOutput(), outPrefix + "residual" + QI::OutExt());
            if (resids) {
                QI::WriteVectorImage(fit_filter->GetResidualsOutput(0), outPrefix + "all_residuals" + QI::OutExt());
            }
            if (its) {
                QI::WriteImage(fit_filter->GetFlagOutput(), outPrefix + "iterations" + QI::OutExt());
            }
            QI_LOG(verbose, "Finished." );
        }
    };
    switch (modelarg.Get()) {
    case 2: {
        QI::TwoPoolModel model{sequences, scale.Get()};
        process(model, "2C_");
    } break;
    case 3: {
        QI::ThreePoolModel model{sequences, scale.Get()};
        process(model, "3C_");
    } break;
    default:
        QI_FAIL("Unknow model specifier: " << modelarg.Get());
    }
    return EXIT_SUCCESS;
}
