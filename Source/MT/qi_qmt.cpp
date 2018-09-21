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

// #define QI_DEBUG_BUILD 1
#include "Macro.h"
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

struct RamaniInnerModel {
    using SequenceType = QI::MTSatSequence;
    using DataType = double;
    using ParameterType = double;
    
    static const int NV = 5;
    static const int NF = 2;
    using VaryingArray = QI_ARRAYN(ParameterType, NV);
    using FixedArray   = QI_ARRAYN(ParameterType, NF);
    static std::array<const std::string, NV> varying_names;
    static std::array<const std::string, NF> fixed_names;
    static const FixedArray fixed_defaults;

    VaryingArray start;
    VaryingArray bounds_lo = VaryingArray::Constant(1.0e-12);
    VaryingArray bounds_hi = VaryingArray::Constant(std::numeric_limits<ParameterType>::infinity());

    QI::Lineshapes lineshape = QI::Lineshapes::Gaussian;
    std::shared_ptr<QI::InterpLineshape> interp = nullptr;

    template<typename Derived>
    auto signal(const Eigen::ArrayBase<Derived> &v,
                            const FixedArray &f,
                            const QI::MTSatSequence *s) const -> QI_ARRAY(typename Derived::Scalar)
    {
        // Use Ramani's notation
        const auto &Rb      = 1.0; // Fix
        const auto &RM0a    = v[0];
        const auto &fterm   = v[1];
        const auto &T2b     = v[2];
        const auto &T1a_T2a = v[3];
        const auto &gM0a    = v[4];
        const auto &f0      = f[0];
        const auto &B1      = f[1];

        QI_ARRAY(typename Derived::Scalar) lsv;
        switch (lineshape) {
        case QI::Lineshapes::Gaussian: lsv = QI::Gaussian((s->sat_f0 + f0), T2b); break;
        case QI::Lineshapes::Lorentzian: lsv = QI::Lorentzian((s->sat_f0 + f0), T2b); break;
        case QI::Lineshapes::SuperLorentzian: lsv = QI::SuperLorentzian((s->sat_f0 + f0), T2b); break;
        case QI::Lineshapes::Interpolated: lsv = (*interp)((s->sat_f0 + f0), T2b); break;
        }

        const auto w_cwpe = (B1 * s->sat_angle / s->pulse.p1) * sqrt(s->pulse.p2 / (s->pulse.Trf * s->TR));
        const auto R_rfb = M_PI * (w_cwpe * w_cwpe) * lsv;
        
        const auto S = gM0a * (Rb*RM0a*fterm + R_rfb + Rb + RM0a) /
            ((RM0a*fterm)*(Rb + R_rfb) + (1.0 + pow(w_cwpe / (2*M_PI*s->sat_f0),2.0)*T1a_T2a)*(R_rfb + Rb + RM0a));
        QI_DBVEC(s->sat_angle);
        QI_DBVEC(s->sat_f0);
        // QI_DBVECT(lsv);
        // QI_DBVECT(w_cwpe);
        QI_DBVECT(v);
        QI_DBVEC(f);
        QI_DBVECT(R_rfb);
        QI_DBVECT(S)
        return S;
    }
};
std::array<const std::string, 5> RamaniInnerModel::varying_names{{"R*M0a"s, "f/(R_a*(1-f))"s, "T2_b"s, "T1_a/T2_a"s, "gM0_a"s}};
std::array<const std::string, 2> RamaniInnerModel::fixed_names{{"f0"s, "B1"s}};
const QI_ARRAYN(double, 2) RamaniInnerModel::fixed_defaults{0.0, 1.0};

struct RamaniOuterModel {
    using SequenceType = QI::MTSatSequence;
    using DataType = double;
    using ParameterType = double;
    
    static const int NV = 6;
    static const int NF = 3;
    using VaryingArray = QI_ARRAYN(ParameterType, NV);
    using FixedArray   = QI_ARRAYN(ParameterType, NF);
    static std::array<const std::string, NV> varying_names;
    static std::array<const std::string, NF> fixed_names;
    static const FixedArray fixed_defaults;

    VaryingArray bounds_lo = VaryingArray::Constant(1.0e-12);
    VaryingArray bounds_hi = VaryingArray::Constant(std::numeric_limits<ParameterType>::infinity());

    QI::Lineshapes lineshape = QI::Lineshapes::Gaussian;
    std::shared_ptr<QI::InterpLineshape> interp = nullptr;

    template<typename Derived>
    auto signal(const Eigen::ArrayBase<Derived> &v,
                            const FixedArray &f,
                            const QI::MTSatSequence *s) const -> QI_ARRAY(typename Derived::Scalar)
    {
        const auto &PD   = v[0];
        const auto &T1_f = v[1];
        const auto &T2_f = v[2];
        const auto &T1_b = 1.0; // Fix
        const auto &T2_b = v[3];
        const auto &k_bf = v[4];
        const auto &f_b  = v[5];
        const auto &f0   = f[0];
        const auto &B1   = f[1];

        QI_DB(PD);
        QI_DB(T1_f);
        QI_DB(T2_f);
        QI_DB(T1_b);
        QI_DB(T2_b);
        QI_DB(k_bf);
        QI_DB(f_b);
        QI_ARRAY(typename Derived::Scalar) lsv;
        switch (lineshape) {
        case QI::Lineshapes::Gaussian: lsv = QI::Gaussian((s->sat_f0 + f0), T2_b); break;
        case QI::Lineshapes::Lorentzian: lsv = QI::Lorentzian((s->sat_f0 + f0), T2_b); break;
        case QI::Lineshapes::SuperLorentzian: lsv = QI::SuperLorentzian((s->sat_f0 + f0), T2_b); break;
        case QI::Lineshapes::Interpolated: lsv = (*interp)((s->sat_f0 + f0), T2_b); break;
        }
        const auto w_cwpe = (B1 * s->sat_angle / s->pulse.p1) * sqrt(s->pulse.p2 / (s->pulse.Trf * s->TR));
        const auto R_rfb = (w_cwpe * w_cwpe) * M_PI * lsv;
        const auto R1_f = 1. / T1_f;
        const auto R1_b = 1. / T1_b;
        
        const auto f_f = 1.0 - f_b;
        const auto k_fb = k_bf * (f_b / f_f);
        const auto M0_f = PD * f_f;
        const auto M0_b = PD * f_b;

        const auto S = M0_f * (R1_b * (k_fb * M0_b / R1_f) + R_rfb + R1_b + k_fb*M0_f) /
            ((k_fb * M0_b / R1_f) * (R1_b + R_rfb) + (1.0 + pow(w_cwpe / (2*M_PI*s->sat_f0),2.0) / (R1_f * T2_f))*(R_rfb + R1_b + k_fb*M0_f));
        QI_DBVEC(s->sat_angle);
        QI_DBVEC(s->sat_f0);
        QI_DBVECT(lsv);
        QI_DBVECT(w_cwpe);
        QI_DBVECT(R_rfb);
        QI_DB(R1_f);
        QI_DB(R1_b);
        QI_DB(f_f);
        QI_DB(k_fb);
        QI_DB(M0_f);
        QI_DB(M0_b);
        QI_DBVECT(S)
        return S;
    }
};
std::array<const std::string, 6> RamaniOuterModel::varying_names{{"PD"s, "T1_f"s, "T2_f"s, "T2_b"s, "k_bf"s, "f_b"s}};
std::array<const std::string, 3> RamaniOuterModel::fixed_names{{"f0"s, "B1"s, "T1_app"}};
const QI_ARRAYN(double, 3) RamaniOuterModel::fixed_defaults{0.0, 1.0, 1.0};

struct RamaniFitFunction : QI::FitFunction<RamaniOuterModel> {
    int max_iterations = 50;
    QI::FitReturnType fit(const std::vector<QI_ARRAY(InputType)> &inputs,
                          const Eigen::ArrayXd &fixed, QI_ARRAYN(OutputType, RamaniOuterModel::NV) &p,
                          ResidualType &residual,
                          std::vector<QI_ARRAY(ResidualType)> &residuals,
                          FlagType &iterations) const override
    {
        const double &scale = inputs[0].maxCoeff();
        if (scale < std::numeric_limits<double>::epsilon()) {
            p = RamaniOuterModel::VaryingArray::Zero();
            residual = 0;
            return std::make_tuple(false, "Maximum data value was not positive");
        }
        const Eigen::ArrayXd data = inputs[0] / scale;

        RamaniInnerModel imodel;
        imodel.lineshape = this->model.lineshape;
        if (imodel.lineshape == QI::Lineshapes::Interpolated) {
            imodel.interp = this->model.interp;
        }
        RamaniInnerModel::VaryingArray inner_p;
        //{"R*M0a"s, "f/(R_a*(1-f))"s, "T2_b"s, "T1_a/T2_a"s, "gM0_a"s}
        inner_p << 10.0, 0.1, 10.e-6, 15., 1.;
        QI_DBVEC(data);
        QI_DBVEC(inner_p);
        ceres::Problem problem;
        using Cost     = QI::ModelCost<RamaniInnerModel>;
        using AutoCost = ceres::AutoDiffCostFunction<Cost, ceres::DYNAMIC, RamaniInnerModel::NV>;
        auto *cost = new Cost(imodel, this->sequence, fixed.head(2), data);
        auto *auto_cost = new AutoCost(cost, this->sequence->size());
        problem.AddResidualBlock(auto_cost, NULL, inner_p.data());
        problem.SetParameterLowerBound(inner_p.data(), 4, this->model.bounds_lo[0] / scale);
        problem.SetParameterUpperBound(inner_p.data(), 4, this->model.bounds_hi[0] / scale);
        for (int i = 0; i < 4; i++) {
            problem.SetParameterLowerBound(inner_p.data(), i, this->model.bounds_lo[i]);
            problem.SetParameterUpperBound(inner_p.data(), i, this->model.bounds_hi[i]);
        }
        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
        options.max_num_iterations = this->max_iterations;
        options.function_tolerance = 1e-5;
        options.gradient_tolerance = 1e-6;
        options.parameter_tolerance = 1e-4;
        options.logging_type = ceres::SILENT;
        ceres::Solve(options, &problem, &summary);
        if (!summary.IsSolutionUsable()) {
            return std::make_tuple(false, summary.FullReport());
        }
        // Convert from the fitted parameters to useful ones
        const auto T1_obs = fixed[2];
        const auto R_obs = 1 / T1_obs;
        const auto &Rb      = 1.0; // Fix
        const auto &RM0a    = inner_p[0];
        const auto &fterm   = inner_p[1];
        const auto &T2b     = inner_p[2];
        const auto &T1a_T2a = inner_p[3];
        const auto &gM0a    = inner_p[4];
        const auto Ra = R_obs / (1.0 + ((RM0a*fterm * (Rb - R_obs))/(Rb - R_obs + RM0a)));
        const auto f = fterm*Ra / (1.0 + fterm*Ra);
        const auto k_bf = RM0a * f / (1.0 - f);
        //{"R*M0a"s, "f/(R_a*(1-f))"s, "T2_b"s, "T1_a/T2_a"s, "gM0_a"s}
        //{"PD"s, "T1_f"s, "T2_f"s, "T2_b"s, "k_bf"s, "f_b"s}
        p[0] = gM0a * scale;
        p[1] = 1.0 / Ra;
        p[2] = p[1] / T1a_T2a;
        p[3] = T2b;
        p[4] = k_bf;
        p[5] = f;
        QI_DBVEC(inner_p);
        QI_DBVEC(p);
        iterations = summary.iterations.size();
        residual = summary.final_cost * scale;
         if (residuals.size() > 0) {
            std::vector<double> r_temp(data.size());
            problem.Evaluate(ceres::Problem::EvaluateOptions(), NULL, &r_temp, NULL, NULL);
            for (size_t i = 0; i < r_temp.size(); i++)
                residuals[0][i] = r_temp[i] * scale;
        }
        return std::make_tuple(true, "");
    }
};



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
    args::ValueFlag<std::string> T1(parser, "T1", "T1 map (seconds) file", {'t', "T1"});
    args::ValueFlag<std::string> lineshape_arg(parser, "LINESHAPE", "Either Gaussian, Lorentzian, Superlorentzian, or a .json file generated by qi_lineshape", {'l', "lineshape"}, "Gaussian");
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
    QI::Lineshapes lineshape;
    std::shared_ptr<QI::InterpLineshape> interp = nullptr;
    if (lineshape_arg.Get() == "Gaussian") {
        QI_LOG(verbose, "Using a Gaussian lineshape");
        lineshape = QI::Lineshapes::Gaussian;
    } else if (lineshape_arg.Get() == "Lorentzian") {
        QI_LOG(verbose, "Using a Lorentzian lineshape");
        lineshape = QI::Lineshapes::Lorentzian;
    } else if (lineshape_arg.Get() == "Superlorentzian") {
        QI_LOG(verbose, "Using a Super-Lorentzian lineshape");
        lineshape = QI::Lineshapes::SuperLorentzian;
    } else {
        QI_LOG(verbose, "Reading lineshape file: " << lineshape_arg.Get());
        rapidjson::Document ls_file = QI::ReadJSON(lineshape_arg.Get());
        interp = std::make_shared<QI::InterpLineshape>(QI::GetMember(ls_file, "lineshape"));
        lineshape = QI::Lineshapes::Interpolated;
    }

    if (simulate) {
        RamaniOuterModel model;
        QI::SimulateModel<RamaniOuterModel, false>(input, model, {&mtsat_sequence}, {f0.Get(), B1.Get()}, {mtsat_path.Get()}, verbose, simulate.Get());
    } else {
        RamaniFitFunction fit;
        fit.model.lineshape = lineshape;
        fit.model.interp = interp;
        fit.sequence = &mtsat_sequence;
        auto fit_filter = itk::ModelFitFilter<RamaniFitFunction>::New();
        fit_filter->SetVerbose(verbose);
        fit_filter->SetFitFunction(&fit);
        fit_filter->SetOutputAllResiduals(resids);
        fit_filter->SetInput(0, QI::ReadVectorImage(mtsat_path.Get(), verbose));
        if (f0) fit_filter->SetFixed(0, QI::ReadImage(f0.Get(), verbose));
        if (B1) fit_filter->SetFixed(1, QI::ReadImage(B1.Get(), verbose));
        if (T1) fit_filter->SetFixed(2, QI::ReadImage(T1.Get(), verbose));
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
        for (int i = 0; i < RamaniOuterModel::NV; i++) {
            QI::WriteImage(fit_filter->GetOutput(i), outPrefix + RamaniOuterModel::varying_names.at(i) + QI::OutExt());
        }
        QI::WriteImage(fit_filter->GetResidualOutput(), outPrefix + "residual" + QI::OutExt());
        if (resids) {
            QI::WriteVectorImage(fit_filter->GetResidualsOutput(0), outPrefix + "all_residuals" + QI::OutExt());
        }
        QI_LOG(verbose, "Finished." );
    }
    return EXIT_SUCCESS;
}
