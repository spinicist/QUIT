/*
 *  qi_esmt.cpp
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "args.hxx"
#include "ceres/ceres.h"

#include "QI/Util.h"
#include "QI/IO.h"
#include "QI/Args.h"
#include "Filters/ApplyAlgorithmFilter.h"

struct EMTCost {
public:
    const Eigen::ArrayXd &G;
    const Eigen::ArrayXd &b;
    const Eigen::ArrayXd &flip;
    const Eigen::ArrayXd &int_omega2;
    const Eigen::ArrayXd &TR;
    const Eigen::ArrayXd &Trf;
    const double T2r;
    const double T2f;
    const double f0_Hz;
    const bool debug;

    template<typename T>
    bool operator() (T const* const* p, T* resids) const {
        typedef Eigen::Array<T, Eigen::Dynamic, 1> ArrayXT;

        const T &M0  = p[0][0];
        const T &F   = p[0][1];
        const T &kf  = p[0][2];
        const T &T1f = p[0][3];
        const T &T1r = T1f; //p[0][3];
        //const double &T2r = 12.0e-6; //p[0][4];

        const ArrayXT E1f = (-TR/T1f).exp();
        const Eigen::ArrayXd E2f = (-TR/T2f).exp();
        const Eigen::ArrayXd E2f_echo = (-TR/(2.0*T2f)).exp();
        const T kr = (F > 0.0) ? (kf / F) : T(0.0);
        //const double T1r = 1.0; // Fixed for now
        const ArrayXT E1r = (-TR/T1r).exp();
        const ArrayXT fk = (-TR*(kf + kr)).exp();

        //const double T2r = 25.e-6; // Fixed for now;
        const double G_gauss = (T2r / sqrt(2.*M_PI))*exp(-pow(2.*M_PI*f0_Hz*T2r,2) / 2.0);
        const Eigen::ArrayXd WT = M_PI * int_omega2 * G_gauss; // # Product of W and Trf to save a division and multiplication
        const Eigen::ArrayXd fw = (-WT).exp();
        const ArrayXT A = 1.0 + F - fw*E1r*(F+fk);
        const ArrayXT B = 1.0 + fk*(F-fw*E1r*(F+1.0));
        const ArrayXT C = F*(1.0-E1r)*(1.0-fk);

        const ArrayXT denom = (A - B*E1f*cos(flip) - (E2f*E2f)*(B*E1f-A*cos(flip)));
        const ArrayXT Gp = M0*E2f_echo*(sin(flip)*((1.0-E1f)*B+C))/denom;
        const ArrayXT bp = (E2f*(A-B*E1f)*(1.0+cos(flip)))/denom;

        Eigen::Map<ArrayXT> r(resids, G.size() + b.size());
        r.head(G.size()) = (G - Gp);
        r.tail(b.size()) = (b - bp);
        if (debug) {
            std::cerr << "M0=" << M0 << "\nF=" << F << "\nkf=" << kf << "\nT1f=" << T1f << "\nT2f=" << T2f << "\nf0_Hz=" << f0_Hz << "\n"
                      << "flip:" << flip.transpose() << "\n"
                      << "int: " << int_omega2.transpose() << "\n"
                      << "TR:  " << TR.transpose() << "\n"
                      << "Trf: " << Trf.transpose() << "\n"
                      << "E1r: " << E1r.transpose() << "\n"
                      << "fk:  " << fk.transpose() << "\n"
                      << "G_gauss: " << G_gauss << "\n"
                      << "WT:  " << WT.transpose() << "\n"
                      << "fw:  " << fw.transpose() << "\n"
                      << "A:   " << A.transpose() << "\n"
                      << "B:   " << B.transpose() << "\n"
                      << "C:   " << C.transpose() << "\n"
                      << "Gp:  " << Gp.transpose() << "\n"
                      << "bp:  " << bp.transpose() << "\n"
                      << "G:   " << G.transpose() << "\n"
                      << "b:   " << b.transpose() << "\n"
                      << "r:   " << r.transpose() << std::endl;
        }
        return true;
    }
};


class EMT : public QI::ApplyF::Algorithm {
public:
    const static size_t NumOutputs = 5;
protected:
    const Eigen::ArrayXd &flips;
    const Eigen::ArrayXd &intB1;
    const Eigen::ArrayXd &TRs;
    const Eigen::ArrayXd &TRFs;
    const double T2r;
    const bool debug;
public:

    EMT(const Eigen::ArrayXd &f, const Eigen::ArrayXd &iB, const Eigen::ArrayXd &tr, const Eigen::ArrayXd &trf, const double T2, const bool d) :
        flips(f), intB1(iB), TRs(tr), TRFs(trf), T2r(T2), debug(d)
    {
    }

    size_t numInputs() const override { return 3; }
    size_t numConsts() const override { return 2; }
    size_t numOutputs() const override { return NumOutputs; }
    size_t dataSize() const override { return (flips.size() * 3); }
    size_t outputSize(const int i) const override { return 1; }
    virtual std::vector<float> defaultConsts() const override {
        std::vector<float> def(2); // B1, f0
        def[0] = 1.0; def[1] = 0.0;
        return def;
    }
    virtual const TOutput &zero(const size_t i) const override { static const float zero = 0.f; return zero; }
    const std::vector<std::string> & names() const {
        static std::vector<std::string> _names = {"M0", "F", "kf", "T1f", "T2f"};
        return _names;
    }
    virtual bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                       std::vector<TOutput> &outputs, TConst &residual,
                       TInput &resids, TIters &its) const override
    {
        const double B1 = consts[0];
        const double f0_Hz = consts[1];

        Eigen::Map<const Eigen::ArrayXf> in_G(inputs[0].GetDataPointer(), inputs[0].Size());
        Eigen::Map<const Eigen::ArrayXf> in_a(inputs[1].GetDataPointer(), inputs[1].Size());
        Eigen::Map<const Eigen::ArrayXf> in_b(inputs[2].GetDataPointer(), inputs[2].Size());

        const double scale = in_G.mean();
        Eigen::ArrayXd G = in_G.cast<double>() / scale;
        Eigen::ArrayXd a = in_a.cast<double>();
        Eigen::ArrayXd b = in_b.cast<double>();
        Eigen::ArrayXd T2fs = (-TRs.cast<double>() / a.log());
        const double T2f = T2fs.mean(); // Different TRs so have to average afterwards

        Eigen::Array<double, 6, 1> p; p << 15.0, 0.05, 5.0, 1.0;
        ceres::Problem problem;
        auto *cost = new ceres::DynamicAutoDiffCostFunction<EMTCost>(new EMTCost{G, b, flips*B1, intB1*B1*B1, TRs, TRFs, T2r, T2f, f0_Hz, debug});
        cost->AddParameterBlock(4);
        cost->SetNumResiduals(G.size() + b.size());
        problem.AddResidualBlock(cost, NULL, p.data());
        const double not_zero = std::nextafter(0.0, 1.0);
        const double not_one  = std::nextafter(1.0, 0.0);
        problem.SetParameterLowerBound(p.data(), 0, 0.1);
        problem.SetParameterUpperBound(p.data(), 0, 20.0);
        problem.SetParameterLowerBound(p.data(), 1, 1e-6);
        problem.SetParameterUpperBound(p.data(), 1, 0.2 - 1e-6);
        problem.SetParameterLowerBound(p.data(), 2, 0.1);
        problem.SetParameterUpperBound(p.data(), 2, 10.0);
        problem.SetParameterLowerBound(p.data(), 3, 0.5);
        problem.SetParameterUpperBound(p.data(), 3, 5.0);
        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
        options.max_num_iterations = 100;
        options.function_tolerance = 1e-7;
        options.gradient_tolerance = 1e-8;
        options.parameter_tolerance = 1e-6;
        if (!debug) options.logging_type = ceres::SILENT;
        ceres::Solve(options, &problem, &summary);
        if (!summary.IsSolutionUsable()) {
            std::cerr << summary.FullReport() << std::endl;
            std::cerr << "Parameters: " << p.transpose() << " T2f: " << T2f << " B1: " << B1 << std::endl;
            std::cerr << "G: " << G.transpose() << std::endl;
            std::cerr << "a: " << a.transpose() << std::endl;
            std::cerr << "b: " << b.transpose() << std::endl;
            return false;
        } else if (debug) {
            std::cout << summary.FullReport() << std::endl;
        }
        outputs[0] = p[0] * scale;
        outputs[1] = p[1];
        outputs[2] = p[2];
        outputs[3] = p[3];
        outputs[4] = T2f;

        residual = summary.final_cost * scale;

        if (resids.Size() > 0) {
            assert(resids.Size() == data.size());
            std::vector<double> r_temp(G.size() + a.size() + b.size());
            problem.Evaluate(ceres::Problem::EvaluateOptions(), NULL, &r_temp, NULL, NULL);
            for (int i = 0; i < G.size(); i++)
                resids[i] = r_temp[i];
            Eigen::ArrayXd as = (-TRs.cast<double>() / T2f).exp();
            for (int i = 0; i < a.size(); i++) {
                resids[i + G.size()] = as[i] - a[i];
            }
            for (int i = 0; i < b.size(); i++) {
                resids[i + G.size() + a.size()] = r_temp[i + G.size()];
            }
        }

        return true;
    }
};

int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates qMT parameters from ellipse parameters.\nInputs are G, a, b.\nhttp://github.com/spinicist/QUIT");
    
    args::Positional<std::string> G_path(parser, "G_FILE", "Input G file");
    args::Positional<std::string> a_path(parser, "a_FILE", "Input a file");
    args::Positional<std::string> b_path(parser, "b_FILE", "Input b file");
    
    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag     debug(parser, "DEBUG", "Output debugging messages", {'d', "debug"});
    args::Flag     noprompt(parser, "NOPROMPT", "Suppress input prompts", {'n', "no-prompt"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "PREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio)", {'b', "B1"});
    args::ValueFlag<std::string> f0(parser, "f0", "f0 map (in Hertz)", {'f', "f0"});
    args::ValueFlag<double> T2r_us(parser, "T2r", "T2r (in microseconds, default 12)", {"T2r"}, 12);
    args::ValueFlag<std::string> subregion(parser, "REGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    args::Flag     all_residuals(parser, "RESIDUALS", "Write out all residuals", {'r',"all_resids"});
    QI::ParseArgs(parser, argc, argv);
    bool prompt = !noprompt;
    
    if (verbose) std::cout << "Opening file: " << QI::CheckPos(G_path) << std::endl;
    auto G = QI::ReadVectorImage<float>(QI::CheckPos(G_path));
    if (verbose) std::cout << "Opening file: " << QI::CheckPos(a_path) << std::endl;
    auto a = QI::ReadVectorImage<float>(QI::CheckPos(a_path));
    if (verbose) std::cout << "Opening file: " << QI::CheckPos(b_path) << std::endl;
    auto b = QI::ReadVectorImage<float>(QI::CheckPos(b_path));

    if (prompt) std::cout << "Enter flip-angles (degrees): ";
    Eigen::ArrayXd flips; QI::ReadArray(std::cin, flips); flips *= M_PI/180.;
    if (prompt) std::cout << "Enter integral B1^2: ";
    Eigen::ArrayXd intB1; QI::ReadArray(std::cin, intB1);
    if (prompt) std::cout << "Enter TRs (seconds): ";
    Eigen::ArrayXd TRs; QI::ReadArray(std::cin, TRs);
    if (prompt) std::cout << "Enter TRFs (seconds): ";
    Eigen::ArrayXd TRFs; QI::ReadArray(std::cin, TRFs);
    if (verbose) {
        std::cout << "T2r " << T2r_us.Get() << "us" << std::endl;
    }
    auto algo = std::make_shared<EMT>(flips, intB1, TRs, TRFs, T2r_us.Get() * 1e-6, debug);

    auto apply = QI::ApplyF::New();
    apply->SetAlgorithm(algo);
    apply->SetPoolsize(threads.Get());
    apply->SetSplitsPerThread(threads.Get()*2);
    apply->SetOutputAllResiduals(all_residuals);
    apply->SetInput(0, G);
    apply->SetInput(1, a);
    apply->SetInput(2, b);
    if (B1) apply->SetConst(0, QI::ReadImage(B1.Get()));
    if (f0) apply->SetConst(1, QI::ReadImage(f0.Get()));
    if (mask) apply->SetMask(QI::ReadImage(mask.Get()));

    apply->SetVerbose(verbose);
    if (subregion) {
        apply->SetSubregion(QI::RegionArg(args::get(subregion)));
    }
    if (verbose) {
        std::cout << "Flips: " << flips.transpose() << std::endl;
        std::cout << "Int B1^2: " << intB1.transpose() << std::endl;
        std::cout << "TR: " << TRs.transpose() << std::endl;
        std::cout << "Trf: " << TRFs.transpose() << std::endl;
        std::cout << "Processing" << std::endl;
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    if (verbose) {
        std::cout << "Elapsed time was " << apply->GetTotalTime() << "s" << std::endl;
        std::cout << "Writing results files." << std::endl;
    }
    std::string outPrefix = outarg.Get() + "EMT_";
    for (int i = 0; i < algo->numOutputs(); i++) {
        std::string outName = outPrefix + algo->names().at(i) + QI::OutExt();
        if (verbose) std::cout << "Writing: " << outName << std::endl;
        QI::WriteImage(apply->GetOutput(i), outName);
    }
    if (verbose) std::cout << "Writing total residual." << std::endl;
    QI::WriteImage(apply->GetResidualOutput(), outPrefix + "residual" + QI::OutExt());
    if (all_residuals) {
        if (verbose) std::cout << "Writing individual residuals." << std::endl;
        QI::WriteScaledVectorImage(apply->GetAllResidualsOutput(), apply->GetOutput(0), outPrefix + "all_residuals" + QI::OutExt());
    }

    if (verbose) std::cout << "Finished." << std::endl;
    return EXIT_SUCCESS;
}
