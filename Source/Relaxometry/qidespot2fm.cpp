/*
 *  despot2fm.cpp
 *
 *  Created by Tobias Wood on 2015/06/03.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <Eigen/Dense>
#include "ceres/ceres.h"

#include "Util.h"
#include "ImageIO.h"
#include "Args.h"
#include "Models.h"
#include "ApplyTypes.h"
#include "SSFPSequence.h"
#include "SequenceCereal.h"

class FMCost : public ceres::CostFunction {
private:
    const Eigen::ArrayXd &m_data;
    double m_T1, m_B1;
    QI::SSFPSequence m_sequence;

public:
    FMCost(const Eigen::ArrayXd &d, const QI::SSFPSequence &s,
           const double T1, const double B1) :
        m_data(d), m_T1(T1), m_B1(B1), m_sequence(s)
    {
        mutable_parameter_block_sizes()->push_back(3);
        set_num_residuals(d.size());
    }

    Eigen::ArrayXd residuals(const Eigen::VectorXd &p) const {
        Eigen::ArrayXd s = QI::One_SSFP_Echo_Magnitude(m_sequence.FA, m_sequence.PhaseInc, m_sequence.TR, p[0], m_T1, p[1], p[2], m_B1);
        Eigen::ArrayXd diff = s - m_data;
        return diff;
    }

    bool Evaluate(double const* const* parameters,
                  double* resids,
                  double** jacobians) const override
    {
        Eigen::Map<const Eigen::Array3d> p(parameters[0]);
        Eigen::Map<Eigen::ArrayXd> r(resids, m_data.size());
        r = residuals(p);
        if (jacobians && jacobians[0]) {
            Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>> j(jacobians[0], m_data.size(), p.size());
            j = QI::One_SSFP_Echo_Derivs(m_sequence.FA, m_sequence.PhaseInc, m_sequence.TR, p[0], m_T1, p[1], p[2], m_B1);
        }
        return true;
    }

};

class LM_FM : public QI::ApplyF::Algorithm {
protected:
    QI::SSFPSequence m_sequence;
    bool m_asymmetric = false, m_debug = false;
public:
    LM_FM(QI::SSFPSequence s, const bool a, const bool d) :
        m_sequence(s), m_asymmetric(a), m_debug(d)
    {}

    size_t numInputs() const override  { return m_sequence.count(); }
    size_t numConsts() const override  { return 2; }
    size_t numOutputs() const override { return 3; }
    size_t dataSize() const override   { return m_sequence.size(); }
    float zero() const override { return 0.f; }
    std::vector<float> defaultConsts() const override {
        std::vector<float> def(2, 1.0f); // T1 & B1
        return def;
    }

    bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
               const TIndex &, // Unused
               std::vector<TOutput> &outputs, TConst &residual,
               TInput &resids, TIterations &its) const override
    {
        const double T1 = consts[0];
        const double B1 = consts[1];
        if (std::isfinite(T1) && (T1 > m_sequence.TR)) {
            // Improve scaling by dividing the PD down to something sensible.
            // This gets scaled back up at the end.
            Eigen::Map<const Eigen::ArrayXf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
            Eigen::ArrayXd data = indata.cast<double>() / indata.maxCoeff();

            std::vector<double> f0_starts = {0, 0.4/m_sequence.TR};
            if (this->m_asymmetric) {
                f0_starts.push_back(0.2/m_sequence.TR);
                f0_starts.push_back(-0.2/m_sequence.TR);
                f0_starts.push_back(-0.4/m_sequence.TR);
            }

            double best = std::numeric_limits<double>::infinity();
            Eigen::Array3d p, bestP;
            ceres::Problem problem;
            problem.AddResidualBlock(new FMCost(data, m_sequence, T1, B1), NULL, p.data());
            problem.SetParameterLowerBound(p.data(), 0, 1.);
            problem.SetParameterLowerBound(p.data(), 1, m_sequence.TR);
            problem.SetParameterUpperBound(p.data(), 1, T1);
            problem.SetParameterLowerBound(p.data(), 2, -0.5/m_sequence.TR);
            problem.SetParameterUpperBound(p.data(), 2,  0.5/m_sequence.TR);
            ceres::Solver::Options options;
            ceres::Solver::Summary summary;
            options.max_num_iterations = 75;
            options.function_tolerance = 1e-6;
            options.gradient_tolerance = 1e-7;
            options.parameter_tolerance = 1e-5;
            if (!m_debug) options.logging_type = ceres::SILENT;
            for (const double &f0 : f0_starts) {
                p = {5., std::max(0.1 * T1, 1.5*m_sequence.TR), f0}; // Yarnykh gives T2 = 0.045 * T1 in brain, but best to overestimate for CSF
                ceres::Solve(options, &problem, &summary);
                if (!summary.IsSolutionUsable()) {
                    std::cerr << summary.FullReport() << std::endl;
                    std::cerr << "T1: " << T1 << " B1: " << B1 << std::endl;
                    std::cerr << "Parameters: " << p.transpose() << std::endl;
                    std::cerr << "Data: " << indata.transpose() << std::endl;
                    return false;
                }
                double r = summary.final_cost;
                if (r < best) {
                    best = r;
                    bestP = p;
                    its = summary.iterations.size();
                }
            }
            if (m_debug) std::cout << summary.FullReport() << std::endl;
            outputs[0] = bestP[0] * indata.maxCoeff();
            outputs[1] = bestP[1];
            outputs[2] = bestP[2];
            
            residual = best * indata.maxCoeff();
            if (resids.Size() > 0) {
                assert(resids.Size() == data.size());
                std::vector<double> r_temp(data.size());
                p = bestP; // Make sure the correct parameters are in the block
                problem.Evaluate(ceres::Problem::EvaluateOptions(), NULL, &r_temp, NULL, NULL);
                for (int i = 0; i < r_temp.size(); i++)
                    resids[i] = r_temp[i];
            }
        } else {
            outputs[0] = 0.;
            outputs[1] = 0.;
            outputs[2] = 0.;
            residual = 0;
            resids.Fill(0.);
            its = 0;
        }
        return true;
    }
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates a T2 map from SSFP data and a T1 map.\nhttp://github.com/spinicist/QUIT");
    
    args::Positional<std::string> t1_path(parser, "T1_MAP", "Input T1 map");
    args::Positional<std::string> ssfp_path(parser, "SSFP_FILE", "Input SSFP file");
    
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio) file", {'b', "B1"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::Flag asym(parser, "ASYM", "Fit +/- off-resonance frequency", {'A', "asym"});
    args::ValueFlag<std::string> subregion(parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    args::Flag debug(parser, "DEBUG", "Output debugging messages", {'d', "debug"});
    args::Flag resids(parser, "RESIDS", "Write out residuals for each data-point", {'r', "resids"});
    QI::ParseArgs(parser, argc, argv, verbose);

    if (verbose) std::cout << "Reading T1 Map from: " << QI::CheckPos(t1_path) << std::endl;
    auto T1 = QI::ReadImage(QI::CheckPos(t1_path));
    if (verbose) std::cout << "Opening SSFP file: " << QI::CheckPos(ssfp_path) << std::endl;
    auto ssfpData = QI::ReadVectorImage<float>(QI::CheckPos(ssfp_path));

    auto ssfp_sequence = QI::ReadSequence<QI::SSFPSequence>(std::cin, verbose);
    auto apply = QI::ApplyF::New();
    std::shared_ptr<LM_FM> algo = std::make_shared<LM_FM>(ssfp_sequence, asym, debug);
    apply->SetVerbose(verbose);
    apply->SetAlgorithm(algo);
    apply->SetOutputAllResiduals(resids);
    if (verbose) std::cout << "Using " << threads.Get() << " threads" << std::endl;
    apply->SetPoolsize(threads.Get());
    apply->SetSplitsPerThread(threads.Get()); // Fairly unbalanced algorithm
    apply->SetInput(0, ssfpData);
    apply->SetConst(0, T1);
    if (B1) apply->SetConst(1, QI::ReadImage(B1.Get()));
    if (mask) apply->SetMask(QI::ReadImage(mask.Get()));
    if (subregion) apply->SetSubregion(QI::RegionArg(args::get(subregion)));
    if (verbose) {
        std::cout << "Processing" << std::endl;
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    if (verbose) {
        std::cout << "Elapsed time was " << apply->GetTotalTime() << "s" << std::endl;
        std::cout << "Writing results files." << std::endl;
    }
    std::string outPrefix = args::get(outarg) + "FM_";
    QI::WriteImage(apply->GetOutput(0), outPrefix + "PD" + QI::OutExt());
    QI::WriteImage(apply->GetOutput(1), outPrefix + "T2" + QI::OutExt());
    QI::WriteImage(apply->GetOutput(2), outPrefix + "f0" + QI::OutExt());
    QI::WriteImage(apply->GetIterationsOutput(), outPrefix + "its" + QI::OutExt());
    QI::WriteScaledImage(apply->GetResidualOutput(), apply->GetOutput(0), outPrefix + "residual" + QI::OutExt());
    if (resids) {
        QI::WriteScaledVectorImage(apply->GetAllResidualsOutput(), apply->GetOutput(0), outPrefix + "all_residuals" + QI::OutExt());
    }
    return EXIT_SUCCESS;
}
