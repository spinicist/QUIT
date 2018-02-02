/*
 *  qdespot1.cpp - Part of QUantitative Imaging Tools
 *
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

#include "ApplyTypes.h"
#include "Models.h"
#include "SPGRSequence.h"
#include "SequenceCereal.h"
#include "Util.h"
#include "Args.h"
#include "ImageIO.h"

//******************************************************************************
// Algorithm Subclasses
//******************************************************************************
class D1Algo : public QI::ApplyF::Algorithm {
public:
    static const size_t DefaultIterations = 15;
protected:
    const std::shared_ptr<QI::Model> m_model = std::make_shared<QI::SCD>();
    QI::SPGRSequence m_sequence;
    size_t m_iterations = DefaultIterations;
    double m_loPD = -std::numeric_limits<double>::infinity();
    double m_hiPD = std::numeric_limits<double>::infinity();
    double m_loT1 = -std::numeric_limits<double>::infinity();
    double m_hiT1 = std::numeric_limits<double>::infinity();

public:
    void setIterations(size_t n) { m_iterations = n; }
    size_t getIterations() { return m_iterations; }
    void setSequence(QI::SPGRSequence &s) { m_sequence = s; }
    void setClampT1(double lo, double hi) { m_loT1 = lo; m_hiT1 = hi; }
    void setClampPD(double lo, double hi) { m_loPD = lo; m_hiPD = hi; }
    size_t numInputs() const override { return m_sequence.count(); }
    size_t numConsts() const override { return 1; }
    size_t numOutputs() const override { return 2; }
    size_t dataSize() const override { return m_sequence.size(); }
    float zero() const override { return 0.f; }
    std::vector<float> defaultConsts() const override {
        // B1
        std::vector<float> def(1, 1.0f);
        return def;
    }
};

class D1LLS : public D1Algo {
public:
    bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
               const TIndex &, // Unused
               std::vector<TOutput> &outputs, TConst &residual,
               TInput &resids, TIterations &its) const override
    {
        Eigen::Map<const Eigen::ArrayXf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
        Eigen::ArrayXd data = indata.cast<double>();
        double B1 = consts[0];
        Eigen::ArrayXd flip = m_sequence.FA * B1;
        Eigen::VectorXd Y = data / flip.sin();
        Eigen::MatrixXd X(Y.rows(), 2);
        X.col(0) = data / flip.tan();
        X.col(1).setOnes();
        Eigen::VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
        outputs[0] = QI::Clamp(b[1] / (1. - b[0]), m_loPD, m_hiPD);
        outputs[1] = QI::Clamp(-m_sequence.TR / log(b[0]), m_loT1, m_hiT1);
        Eigen::ArrayXd theory = QI::One_SPGR(m_sequence.FA, m_sequence.TR, outputs[0], outputs[1], B1).array().abs();
        Eigen::ArrayXf r = (data.array() - theory).cast<float>();
        residual = sqrt(r.square().sum() / r.rows());
        resids = itk::VariableLengthVector<float>(r.data(), r.rows());
        its = 1;
        return true;
    }
};

class D1WLLS : public D1Algo {
public:
    bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
               const TIndex &, // Unused
               std::vector<TOutput> &outputs, TConst &residual,
               TInput &resids, TIterations &its) const override
    {
        Eigen::Map<const Eigen::ArrayXf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
        Eigen::ArrayXd data = indata.cast<double>();
        double B1 = consts[0];
        Eigen::ArrayXd flip = m_sequence.FA * B1;
        Eigen::VectorXd Y = data / flip.sin();
        Eigen::MatrixXd X(Y.rows(), 2);
        X.col(0) = data / flip.tan();
        X.col(1).setOnes();
        Eigen::Vector2d b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
        Eigen::Array2d out;
        out[1] = -m_sequence.TR / log(b[0]);
        out[0] = b[1] / (1. - b[0]);
        for (its = 0; its < m_iterations; its++) {
            Eigen::VectorXd W = (flip.sin() / (1. - (exp(-m_sequence.TR/outputs[1])*flip.cos()))).square();
            b = (X.transpose() * W.asDiagonal() * X).partialPivLu().solve(X.transpose() * W.asDiagonal() * Y);
            Eigen::Array2d newOut;
            newOut[1] = -m_sequence.TR / log(b[0]);
            newOut[0] = b[1] / (1. - b[0]);
            if (newOut.isApprox(out))
                break;
            else
                out = newOut;
        }
        outputs[0] = QI::Clamp(out[0], m_loPD, m_hiPD);
        outputs[1] = QI::Clamp(out[1], m_loT1, m_hiT1);
        Eigen::ArrayXd theory = QI::One_SPGR(m_sequence.FA, m_sequence.TR, outputs[0], outputs[1], B1).array().abs();
        Eigen::ArrayXf r = (data.array() - theory).cast<float>();
        residual = sqrt(r.square().sum() / r.rows());
        resids = itk::VariableLengthVector<float>(r.data(), r.rows());
        return true;
    }
};

class T1Cost : public ceres::CostFunction {
protected:
    const QI::SPGRSequence m_seq;
    const Eigen::ArrayXd m_data;
    const double m_B1;

public:
    T1Cost(const QI::SPGRSequence cs, const Eigen::ArrayXd &data, const double B1) :
        m_seq(cs), m_data(data), m_B1(B1)
    {
        mutable_parameter_block_sizes()->push_back(2);
        set_num_residuals(data.size());
    }

    bool Evaluate(double const* const* parameters,
                  double* resids,
                  double** jacobians) const override
    {
        Eigen::Map<const Eigen::Array2d> p(parameters[0]);
        Eigen::Map<Eigen::ArrayXd> r(resids, m_data.size());
        Eigen::ArrayXd s = QI::One_SPGR(m_seq.FA, m_seq.TR, p[0], p[1], m_B1).array().abs();
        r = s - m_data;
        if (jacobians && jacobians[0]) {
            Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>> j(jacobians[0], m_data.size(), p.size());
            j = QI::One_SPGR_Magnitude_Derivs(m_seq.FA, m_seq.TR, p[0], p[1], m_B1);
        }
        return true;
    }


};

class D1NLLS : public D1Algo {
public:
    D1NLLS() {
        m_loT1 = 1e-6;
        m_loPD = 1e-6;
    }

    bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
               const TIndex &, // Unused
               std::vector<TOutput> &outputs, TConst &residual,
               TInput &resids, TIterations &its) const override
    {
        Eigen::Map<const Eigen::ArrayXf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
        const double B1 = consts[0];
        const double scale = indata.maxCoeff();
        if (scale < 0) {
            outputs[0] = 0;
            outputs[1] = 0;
            residual = 0;
            return false;
        }
        const Eigen::ArrayXd data = indata.cast<double>() / scale;
        Eigen::Array2d p; p << 10., 1.;
        ceres::Problem problem;
        problem.AddResidualBlock(new T1Cost(m_sequence, data, B1), NULL, p.data());
        problem.SetParameterLowerBound(p.data(), 0, m_loPD / scale);
        problem.SetParameterUpperBound(p.data(), 0, m_hiPD / scale);
        problem.SetParameterLowerBound(p.data(), 1, m_loT1);
        problem.SetParameterUpperBound(p.data(), 1, m_hiT1);
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
        
        outputs[0] = p[0] * indata.maxCoeff();
        outputs[1] = p[1];
        if (!summary.IsSolutionUsable()) {
            std::cout << summary.FullReport() << std::endl;
        }
        its = summary.iterations.size();
        residual = summary.final_cost * indata.maxCoeff();
        if (resids.Size() > 0) {
            assert(resids.Size() == data.size());
            std::vector<double> r_temp(data.size());
            problem.Evaluate(ceres::Problem::EvaluateOptions(), NULL, &r_temp, NULL, NULL);
            for (int i = 0; i < r_temp.size(); i++)
                resids[i] = r_temp[i];
        }
        return true;
    }
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates T1 maps from SPGR data\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> spgr_path(parser, "SPGR FILE", "Path to SPGR data");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio) file", {'b', "B1"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> subregion(parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    args::Flag resids(parser, "RESIDS", "Write out residuals for each data-point", {'r', "resids"});
    args::ValueFlag<char> algorithm(parser, "ALGO", "Choose algorithm (l/w/n)", {'a',"algo"}, 'l');
    args::ValueFlag<int> its(parser, "ITERS", "Max iterations for WLLS/NLLS (default 15)", {'i',"its"}, 15);
    args::ValueFlag<float> clampPD(parser, "CLAMP PD", "Clamp PD between 0 and value", {'p',"clampPD"}, std::numeric_limits<float>::infinity());
    args::ValueFlag<float> clampT1(parser, "CLAMP T1", "Clamp T1 between 0 and value", {'t',"clampT2"}, std::numeric_limits<float>::infinity());
    QI::ParseArgs(parser, argc, argv);

    if (verbose) std::cout << "Opening SPGR file: " << QI::CheckPos(spgr_path) << std::endl;
    auto data = QI::ReadVectorImage<float>(QI::CheckPos(spgr_path));
    std::shared_ptr<D1Algo> algo;
    switch (algorithm.Get()) {
        case 'l': algo = std::make_shared<D1LLS>();  if (verbose) std::cout << "LLS algorithm selected." << std::endl; break;
        case 'w': algo = std::make_shared<D1WLLS>(); if (verbose) std::cout << "WLLS algorithm selected." << std::endl; break;
        case 'n': algo = std::make_shared<D1NLLS>(); if (verbose) std::cout << "NLLS algorithm selected." << std::endl; break;
    }
    algo->setIterations(its.Get());
    if (clampPD) algo->setClampPD(1e-6, clampPD.Get());
    if (clampT1) algo->setClampT1(1e-6, clampT1.Get());
    auto spgrSequence = QI::ReadSequence<QI::SPGRSequence>(std::cin, verbose);
    algo->setSequence(spgrSequence);
    auto apply = QI::ApplyF::New();
    apply->SetVerbose(verbose);
    apply->SetAlgorithm(algo);
    apply->SetOutputAllResiduals(resids);
    apply->SetPoolsize(threads.Get());
    apply->SetSplitsPerThread(threads.Get()); // Unbalanced algorithm
    apply->SetInput(0, data);
    if (B1) apply->SetConst(0, QI::ReadImage(B1.Get()));
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
    std::string outPrefix = outarg.Get() + "D1_";
    QI::WriteImage(apply->GetOutput(0), outPrefix + "PD" + QI::OutExt());
    QI::WriteImage(apply->GetOutput(1), outPrefix + "T1" + QI::OutExt());
    QI::WriteScaledImage(apply->GetResidualOutput(), apply->GetOutput(0), outPrefix + "residual" + QI::OutExt());
    if (resids) {
        QI::WriteScaledVectorImage(apply->GetAllResidualsOutput(), apply->GetOutput(0), outPrefix + "all_residuals" + QI::OutExt());
    }
    if (its) {
        QI::WriteImage(apply->GetIterationsOutput(), outPrefix + "iterations" + QI::OutExt());
    }
    if (verbose) std::cout << "Finished." << std::endl;
    return EXIT_SUCCESS;
}
