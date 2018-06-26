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

#include <iostream>
#include <Eigen/Dense>
#include "ceres/ceres.h"

#include "Util.h"
#include "SPGRSequence.h"
#include "MPRAGESequence.h"
#include "SequenceCereal.h"
#include "Args.h"
#include "ImageIO.h"
#include "ApplyTypes.h"

class SPGRCost : public ceres::CostFunction {
protected:
    const QI::SPGRSequence &m_seq;
    const Eigen::ArrayXd m_data;

public:
    SPGRCost(const QI::SPGRSequence &s, const Eigen::ArrayXd &data) :
        m_seq(s), m_data(data)
    {
        mutable_parameter_block_sizes()->push_back(3);
        set_num_residuals(data.size());
    }

    bool Evaluate(double const* const* p,
                  double* resids,
                  double** jacobians) const override
    {
        const double &M0 = p[0][0];
        const double &T1 = p[0][1];
        const double &B1 = p[0][2];

        const Eigen::ArrayXd sa = sin(B1 * m_seq.FA);
        const Eigen::ArrayXd ca = cos(B1 * m_seq.FA);
        const double E1 = exp(-m_seq.TR / T1);
        const Eigen::ArrayXd denom = (1.-E1*ca);
        
        Eigen::Map<Eigen::ArrayXd> r(resids, m_data.size());
        r = M0*sa*(1-E1)/denom - m_data;
        
        // std::cout << "SPGR RESIDS" << std::endl;
        // std::cout << r.transpose() << std::endl;
        if (jacobians && jacobians[0]) {
            Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>> j(jacobians[0], m_data.size(), 3);
            j.col(0) = (1-E1)*sa/denom;
            j.col(1) = E1*M0*m_seq.TR*(ca-1.)*sa/((denom*T1).square());
            j.col(2) = M0*m_seq.FA*(1.-E1)*(ca-E1)/denom.square();
        }
        return true;
    }
};

// Use AutoDiff for this
class IRCostFunction  {
protected:
    const QI::MPRAGESequence &m_seq;
    const Eigen::ArrayXd m_data;

public:
    IRCostFunction(const QI::MPRAGESequence &s, const Eigen::ArrayXd &data) :
        m_seq(s), m_data(data)
    {
    }

    template<typename T> bool operator() (const T* const p1, T* r) const
    {
        const T &M0 = p1[0];
        const T &T1 = p1[1];
        const T &B1 = p1[2];
        const double eta = -1.0; // Inversion efficiency defined as -1 < eta < 0

        const double TIs = m_seq.TI - m_seq.TR*m_seq.k0; // Adjust TI for k0
        const T T1s = 1. / (1./T1 - log(cos(m_seq.FA * B1))/m_seq.TR);
        const T M0s = M0 * (1. - exp(-m_seq.TR/T1)) / (1. - exp(-m_seq.TR/T1s));
        const T A_1 = M0s*(1. - exp(-(m_seq.ETL*m_seq.TR)/T1s));

        const T A_2 = M0*(1. - exp(-m_seq.TD/T1));
        const T A_3 = M0*(1. - exp(-TIs/T1));
        const T B_1 = exp(-(m_seq.ETL*m_seq.TR)/T1s);
        const T B_2 = exp(-m_seq.TD/T1);
        const T B_3 = eta*exp(-TIs/T1);

        const T A = A_3 + A_2*B_3 + A_1*B_2*B_3;
        const T B = B_1*B_2*B_3;
        const T M1 = A / (1. - B);

        r[0] = m_data[0] - (M0s + (M1 - M0s)*exp(-(m_seq.k0*m_seq.TR)/T1s)) * sin(m_seq.FA * B1);
        return true;
    }
};

class HIFIAlgo : public QI::ApplyF::Algorithm {
private:
    const QI::SPGRSequence &m_spgr;
    const QI::MPRAGESequence &m_mprage;
    double m_lo = 0;
    double m_hi = std::numeric_limits<double>::infinity();
public:
    HIFIAlgo(const QI::SPGRSequence &s, const QI::MPRAGESequence &m, const float hi) :
        m_spgr(s), m_mprage(m), m_hi(hi)
    {}
    size_t numInputs() const override  { return 2; }
    size_t numConsts() const override  { return 0; }
    size_t numOutputs() const override { return 3; }
    size_t dataSize() const override   { return m_spgr.size() + m_mprage.size(); }
    float zero() const override { return 0.f; }

    std::vector<float> defaultConsts() const override {
        // No constants for HIFI
        std::vector<float> def(0);
        return def;
    }

    TStatus apply(const std::vector<TInput> &inputs, const std::vector<TConst> &, // No constants, remove name to silence compiler warnings
               const TIndex &, // Unused
               std::vector<TOutput> &outputs, TConst &residual,
               TInput &resids, TIterations &its) const override
    {
        Eigen::Map<const Eigen::ArrayXf> spgr_in(inputs[0].GetDataPointer(), inputs[0].Size());
        Eigen::Map<const Eigen::ArrayXf> ir_in(inputs[1].GetDataPointer(), inputs[1].Size());
        double scale = std::max(spgr_in.maxCoeff(), ir_in.maxCoeff());
        if (scale < std::numeric_limits<double>::epsilon()) {
            outputs[0] = 0;
            outputs[1] = 0;
            outputs[2] = 0;
            residual = 0;
            return std::make_tuple(false, "Maximum data value was zero or less");
        }
        const Eigen::ArrayXd spgr_data = spgr_in.cast<double>() / scale;
        const Eigen::ArrayXd ir_data = ir_in.cast<double>() / scale;
        double spgr_pars[] = {10., 1., 1.}; // PD, T1, B1
        ceres::Problem problem;
        problem.AddResidualBlock(new SPGRCost(m_spgr, spgr_data), NULL, spgr_pars);
        ceres::CostFunction *IRCost = new ceres::AutoDiffCostFunction<IRCostFunction, 1, 3>(new IRCostFunction(m_mprage, ir_data));
        problem.AddResidualBlock(IRCost, NULL, spgr_pars);
        problem.SetParameterLowerBound(spgr_pars, 0, 1.);
        problem.SetParameterLowerBound(spgr_pars, 1, 0.001);
        problem.SetParameterUpperBound(spgr_pars, 1, 5.0);
        problem.SetParameterLowerBound(spgr_pars, 2, 0.1);
        problem.SetParameterUpperBound(spgr_pars, 2, 2.0);
        ceres::Solver::Options options;
        ceres::Solver::Summary summary;
        options.max_num_iterations = 50;
        options.function_tolerance = 1e-5;
        options.gradient_tolerance = 1e-6;
        options.parameter_tolerance = 1e-4;
        // options.check_gradients = true;
        options.logging_type = ceres::SILENT;
        // std::cout << "START P: " << p.transpose() << std::endl;
        ceres::Solve(options, &problem, &summary);
        
        outputs[0] = spgr_pars[0] * scale;
        outputs[1] = QI::Clamp(spgr_pars[1], m_lo, m_hi);
        outputs[2] = spgr_pars[2];
        if (!summary.IsSolutionUsable()) {
            return std::make_tuple(false, summary.FullReport());
        }
        its = summary.iterations.size();
        residual = summary.final_cost * scale;
        if (resids.Size() > 0) {
            std::vector<double> r_temp(spgr_data.size() + 1);
            problem.Evaluate(ceres::Problem::EvaluateOptions(), NULL, &r_temp, NULL, NULL);
            for (size_t i = 0; i < r_temp.size(); i++) {
                resids[i] = r_temp[i];
            }
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
    args::Positional<std::string> ir_path(parser, "IRSPGR_FILE", "Input IR-SPGR or MP-RAGE file");
    
    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag     all_resids(parser, "ALL RESIDUALS", "Output individual residuals in addition to the Sum-of-Squares", {'r',"resids"});
    args::ValueFlag<float> clamp(parser, "CLAMP", "Clamp output T1 values to this value", {'c', "clamp"}, std::numeric_limits<float>::infinity());
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> subregion(parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    args::Flag resids(parser, "RESIDS", "Write out residuals for each data-point", {'r', "resids"});
    QI::ParseArgs(parser, argc, argv, verbose);

    if (verbose) std::cout << "Reading SPGR file: " << QI::CheckPos(spgr_path) << std::endl;
    auto spgrImg = QI::ReadVectorImage(QI::CheckPos(spgr_path));
    cereal::JSONInputArchive input(std::cin);
    auto spgr_sequence = QI::ReadSequence<QI::SPGRSequence>(input, verbose);
    if (verbose) std::cout << "Reading MPRAGE file: " << QI::CheckPos(ir_path) << std::endl;
    auto irImg = QI::ReadVectorImage(QI::CheckPos(ir_path));
    auto ir_sequence = QI::ReadSequence<QI::MPRAGESequence>(input, verbose);

    auto apply = QI::ApplyF::New();
    auto hifi = std::make_shared<HIFIAlgo>(spgr_sequence, ir_sequence, clamp.Get());
    apply->SetAlgorithm(hifi);
    apply->SetOutputAllResiduals(all_resids);
    apply->SetPoolsize(threads.Get());
    apply->SetSplitsPerThread(threads.Get());
    apply->SetVerbose(verbose);
    apply->SetInput(0, spgrImg);
    apply->SetInput(1, irImg);
    if (subregion) apply->SetSubregion(QI::RegionArg(args::get(subregion)));
    if (mask) apply->SetMask(QI::ReadImage(mask.Get()));
    if (verbose) std::cout << "Processing..." << std::endl;
    apply->Update();
    if (verbose) {
        std::cout << "Elapsed time was " << apply->GetTotalTime() << "s" << std::endl;
        std::cout << "Writing results files." << std::endl;
    }
    std::string out_prefix = args::get(outarg) + "HIFI_";
    QI::WriteImage(apply->GetOutput(0), out_prefix + "PD" + QI::OutExt());
    QI::WriteImage(apply->GetOutput(1), out_prefix + "T1" + QI::OutExt());
    QI::WriteImage(apply->GetOutput(2), out_prefix + "B1" + QI::OutExt());
    QI::WriteImage(apply->GetResidualOutput(), out_prefix + "residual"  + QI::OutExt());
    if (all_resids) {
        QI::WriteVectorImage(apply->GetAllResidualsOutput(), out_prefix + "all_residuals" + QI::OutExt());
    }
    if (verbose) std::cout << "Finished." << std::endl;
    return EXIT_SUCCESS;
}
