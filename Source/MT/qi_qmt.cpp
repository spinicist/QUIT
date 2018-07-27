/*
 *  qi_qmt.cpp - Part of QUantitative Imaging Tools
 *
 *  Copyright (c) 2018 Tobias Wood, Samuel Hurley, Erika Raven
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>

#include <Eigen/Core>
#include "ceres/ceres.h"

#include "Model.h"
#include "FitFunction.h"
#include "ModelFitFilter.h"
#include "SimulateModel.h"
#include "MTSatSequence.h"
#include "Lineshape.h"
#include "Util.h"
#include "Args.h"
#include "ImageIO.h"

using namespace std::literals;

struct RamaniModel {
    using SequenceType = QI::MTSatSequence;
    using DataType = double;
    using ParameterType = double;
    
    static const int NV = 7;
    static const int NF = 2;
    using VaryingArray = QI_ARRAYN(ParameterType, NV);
    using FixedArray   = QI_ARRAYN(ParameterType, NF);
    static std::array<const std::string, NV> varying_names;
    static std::array<const std::string, NF> fixed_names;
    static const FixedArray fixed_defaults;

    VaryingArray bounds_lo = VaryingArray::Constant(1.0e-12);
    VaryingArray bounds_hi = VaryingArray::Constant(std::numeric_limits<ParameterType>::infinity());

    QI::Lineshapes::Gaussian ls;

    template<typename Derived>
    auto signal(const Eigen::ArrayBase<Derived> &v,
                            const FixedArray &f,
                            const QI::MTSatSequence *s) const -> QI_ARRAY(typename Derived::Scalar)
    {
        const auto &PD   = v[0];
        const auto &T1_f = v[1];
        const auto &T2_f = v[2];
        const auto &T1_b = v[3];
        const auto &T2_b = v[4];
        const auto &k_bf = v[5];
        const auto &f_b  = v[6];
        const auto &f0   = f[0];
        const auto &B1   = f[1];

        const auto W = (B1*s->sat_angle).square()*ls.value((s->sat_f0 + f0), T2_b);
        const auto R1f = 1. / T1_f;
        const auto R1r = 1. / T1_b;
        
        // F is M0r/M0b
        const auto F = f_b / (1.0 - f_b);
        const auto kr = k_bf/F;
        const auto S = PD * F * ( R1r*kr/R1f + W + R1r + kr ) /
                        ( k_bf*(R1r + W) + (1.0 + ((B1*s->sat_angle)/(2.*M_PI*(s->sat_f0 + f0))).square()*(T1_f/T2_f))*(W+R1r+kr));
        return S;
    }
};
std::array<const std::string, 7> RamaniModel::varying_names{{"PD"s, "T1_f"s, "T2_f"s, "T1_b"s, "T2_b"s, "k_bf"s, "f_b"s}};
std::array<const std::string, 2> RamaniModel::fixed_names{{"f0"s, "B1"s}};
const QI_ARRAYN(double, 2) RamaniModel::fixed_defaults{0.0, 1.0};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates qMT maps from Gradient Echo Saturation data\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> mtsat_path(parser, "MTSAT FILE", "Path to MT-Sat data");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, QI::GetDefaultThreads());
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> f0(parser, "f0", "f0 map (Hz) file", {'f', "f0"});
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio) file", {'b', "B1"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> subregion(parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    args::Flag resids(parser, "RESIDS", "Write out residuals for each data-point", {'r', "resids"});
    args::ValueFlag<char> algorithm(parser, "ALGO", "Choose algorithm (l/w/n)", {'a',"algo"}, 'l');
    args::ValueFlag<std::string> seq_arg(parser, "FILE", "Read JSON input from file instead of stdin", {"file"});
    args::ValueFlag<float> simulate(parser, "SIMULATE", "Simulate sequence instead of fit_filterting model (argument is noise level)", {"simulate"}, 0.0);
    QI::ParseArgs(parser, argc, argv, verbose, threads);
    QI::CheckPos(mtsat_path);
    QI_LOG(verbose, "Reading sequence information");
    rapidjson::Document input = seq_arg ? QI::ReadJSON(seq_arg.Get()) : QI::ReadJSON(std::cin);
    QI::MTSatSequence mtsat_sequence(QI::GetMember(input, "MTSat"));
    if (simulate) {
        RamaniModel model;
        QI::SimulateModel<RamaniModel, false>(input, model, {&mtsat_sequence}, {f0.Get(), B1.Get()}, {mtsat_path.Get()}, verbose, simulate.Get());
    } else {
        QI::ScaledNLLSFitFunction<RamaniModel> fit;
        fit.sequence = &mtsat_sequence;
        auto fit_filter = itk::ModelFitFilter<QI::ScaledNLLSFitFunction<RamaniModel>>::New();
        fit_filter->SetVerbose(verbose);
        fit_filter->SetFitFunction(&fit);
        fit_filter->SetOutputAllResiduals(resids);
        fit_filter->SetInput(0, QI::ReadVectorImage(mtsat_path.Get(), verbose));
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
        std::string outPrefix = outarg.Get() + "QMT_";
        for (int i = 0; i < RamaniModel::NV; i++) {
            QI::WriteImage(fit_filter->GetOutput(i), outPrefix + RamaniModel::varying_names.at(i) + QI::OutExt());
        }
        QI::WriteImage(fit_filter->GetResidualOutput(), outPrefix + "residual" + QI::OutExt());
        if (resids) {
            QI::WriteVectorImage(fit_filter->GetResidualsOutput(0), outPrefix + "all_residuals" + QI::OutExt());
        }
        QI_LOG(verbose, "Finished." );
    }
    return EXIT_SUCCESS;
}
