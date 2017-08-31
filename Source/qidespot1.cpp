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

#include "Filters/ApplyAlgorithmFilter.h"
#include "QI/Models/Models.h"
#include "QI/Sequences/Sequences.h"
#include "QI/Util.h"
#include "QI/Option.h"

using namespace std;
using namespace Eigen;

//******************************************************************************
// Algorithm Subclasses
//******************************************************************************
class D1Algo : public QI::ApplyF::Algorithm {
public:
    static const size_t DefaultIterations = 15;
protected:
    const shared_ptr<QI::Model> m_model = make_shared<QI::SCD>();
    shared_ptr<QI::SPGRSimple> m_sequence;
    size_t m_iterations = DefaultIterations;
    double m_loPD = -numeric_limits<double>::infinity();
    double m_hiPD = numeric_limits<double>::infinity();
    double m_loT1 = -numeric_limits<double>::infinity();
    double m_hiT1 = numeric_limits<double>::infinity();

public:
    void setIterations(size_t n) { m_iterations = n; }
    size_t getIterations() { return m_iterations; }
    void setSequence(shared_ptr<QI::SPGRSimple> &s) { m_sequence = s; }
    void setClampT1(double lo, double hi) { m_loT1 = lo; m_hiT1 = hi; }
    void setClampPD(double lo, double hi) { m_loPD = lo; m_hiPD = hi; }
    size_t numInputs() const override { return m_sequence->count(); }
    size_t numConsts() const override { return 1; }
    size_t numOutputs() const override { return 2; }
    size_t dataSize() const override { return m_sequence->size(); }
    const float &zero(const size_t i) const override { static float zero = 0; return zero; }
    virtual std::vector<float> defaultConsts() const override {
        // B1
        std::vector<float> def(1, 1.0f);
        return def;
    }
};

class D1LLS : public D1Algo {
public:
    virtual bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                  std::vector<TOutput> &outputs, TConst &residual,
                  TInput &resids, TIters &its) const override
    {
        Eigen::Map<const Eigen::ArrayXf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
        ArrayXd data = indata.cast<double>();
        double B1 = consts[0];
        ArrayXd flip = m_sequence->flip() * B1;
        VectorXd Y = data / flip.sin();
        MatrixXd X(Y.rows(), 2);
        X.col(0) = data / flip.tan();
        X.col(1).setOnes();
        VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
        outputs[0] = QI::Clamp(b[1] / (1. - b[0]), m_loPD, m_hiPD);
        outputs[1] = QI::Clamp(-m_sequence->TR() / log(b[0]), m_loT1, m_hiT1);
        ArrayXd theory = QI::One_SPGR(m_sequence->flip(), m_sequence->TR(), outputs[0], outputs[1], B1).array().abs();
        ArrayXf r = (data.array() - theory).cast<float>();
        residual = sqrt(r.square().sum() / r.rows());
        resids = itk::VariableLengthVector<float>(r.data(), r.rows());
        its = 1;
        return true;
    }
};

class D1WLLS : public D1Algo {
public:
    virtual bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                  std::vector<TOutput> &outputs, TConst &residual,
                  TInput &resids, TIters &its) const override
    {
        Eigen::Map<const Eigen::ArrayXf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
        ArrayXd data = indata.cast<double>();
        double B1 = consts[0];
        ArrayXd flip = m_sequence->flip() * B1;
        VectorXd Y = data / flip.sin();
        MatrixXd X(Y.rows(), 2);
        X.col(0) = data / flip.tan();
        X.col(1).setOnes();
        Vector2d b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
        Array2d out;
        out[1] = -m_sequence->TR() / log(b[0]);
        out[0] = b[1] / (1. - b[0]);
        for (its = 0; its < m_iterations; its++) {
            VectorXd W = (flip.sin() / (1. - (exp(-m_sequence->TR()/outputs[1])*flip.cos()))).square();
            b = (X.transpose() * W.asDiagonal() * X).partialPivLu().solve(X.transpose() * W.asDiagonal() * Y);
            Array2d newOut;
            newOut[1] = -m_sequence->TR() / log(b[0]);
            newOut[0] = b[1] / (1. - b[0]);
            if (newOut.isApprox(out))
                break;
            else
                out = newOut;
        }
        outputs[0] = QI::Clamp(out[0], m_loPD, m_hiPD);
        outputs[1] = QI::Clamp(out[1], m_loT1, m_hiT1);
        ArrayXd theory = QI::One_SPGR(m_sequence->flip(), m_sequence->TR(), outputs[0], outputs[1], B1).array().abs();
        ArrayXf r = (data.array() - theory).cast<float>();
        residual = sqrt(r.square().sum() / r.rows());
        resids = itk::VariableLengthVector<float>(r.data(), r.rows());
        return true;
    }
};

class T1Cost : public ceres::CostFunction {
protected:
    const shared_ptr<QI::SequenceBase> m_seq;
    const ArrayXd m_data;
    const double m_B1;

public:
    T1Cost(const shared_ptr<QI::SequenceBase> cs, const ArrayXd &data, const double B1) :
        m_seq(cs), m_data(data), m_B1(B1)
    {
        mutable_parameter_block_sizes()->push_back(2);
        set_num_residuals(data.size());
    }

    virtual bool Evaluate(double const* const* parameters,
                          double* resids,
                          double** jacobians) const override
    {
        Eigen::Map<const Eigen::Array2d> p(parameters[0]);
        Eigen::Map<Eigen::ArrayXd> r(resids, m_data.size());
        ArrayXd s = QI::One_SPGR(m_seq->flip(), m_seq->TR(), p[0], p[1], m_B1).array().abs();
        r = s - m_data;
        // std::cout << "p: " << p.transpose() << std::endl;
        // std::cout << "s: " << s.transpose() << std::endl;
        // std::cout << "d: " << m_data.transpose() << std::endl;
        // std::cout << "r: " << r.transpose() << std::endl;
        if (jacobians && jacobians[0]) {
            Eigen::Map<Eigen::Matrix<double, -1, -1, RowMajor>> j(jacobians[0], m_data.size(), p.size());
            j = QI::One_SPGR_Magnitude_Derivs(m_seq->flip(), m_seq->TR(), p[0], p[1], m_B1);
            // std::cout << "j" << std::endl << j << std::endl;
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

    virtual bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                       std::vector<TOutput> &outputs, TConst &residual,
                       TInput &resids, TIters &its) const override
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
        const ArrayXd data = indata.cast<double>() / scale;
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
            vector<double> r_temp(data.size());
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
    QI::OptionList opts("Usage is: qidespot1 [options] spgr_input");
    QI::Switch all_residuals('r',"resids","Write out per flip-angle residuals", opts);
    QI::Option<int> num_threads(4,'T',"threads","Use N threads (default=4, 0=hardware limit)", opts);
    QI::Option<int> its(15,'i',"its","Max iterations for WLLS/NLLS (default 15)", opts);
    QI::Option<float> clampPD(std::numeric_limits<float>::infinity(),'p',"clampPD","Clamp PD between 1e-6 and value", opts);
    QI::Option<float> clampT1(std::numeric_limits<float>::infinity(),'t',"clampT1","Clamp T1 between 1e-6 and value", opts);
    QI::ImageOption<QI::VolumeF> mask('m', "mask", "Mask input with specified file", opts);
    QI::ImageOption<QI::VolumeF> B1('b', "B1", "B1 Map file (ratio)", opts);
    QI::EnumOption algorithm("lwn",'l','a',"algo","Choose algorithm (l/w/n)", opts);
    QI::Option<std::string> outPrefix("", 'o', "out","Add a prefix to output filenames", opts);
    QI::Switch suppress('n',"no-prompt","Suppress input prompts", opts);
    QI::Switch verbose('v',"verbose","Print more information", opts);
    QI::Help help(opts);
    std::deque<std::string> nonopts = opts.parse(argc, argv);
    if (nonopts.size() != 1) {
        std::cerr << opts << std::endl;
        std::cerr << "No input filename specified." << std::endl;
        return EXIT_FAILURE;
    }

    const std::string &inputFilename = nonopts.front();
    if (*verbose) cout << "Opening SPGR file: " << inputFilename << endl;
    auto data = QI::ReadVectorImage<float>(inputFilename);
    shared_ptr<QI::SPGRSimple> spgrSequence = make_shared<QI::SPGRSimple>(cin, !(*suppress));
    if (*verbose) cout << *spgrSequence;
    shared_ptr<D1Algo> algo;
    switch (*algorithm) {
        case 'l': algo = make_shared<D1LLS>();  if (*verbose) cout << "LLS algorithm selected." << endl; break;
        case 'w': algo = make_shared<D1WLLS>(); if (*verbose) cout << "WLLS algorithm selected." << endl; break;
        case 'n': algo = make_shared<D1NLLS>(); if (*verbose) cout << "NLLS algorithm selected." << endl; break;
    }
    algo->setIterations(*its);
    if (isfinite(*clampPD))
        algo->setClampPD(1e-6, *clampPD);
    if (isfinite(*clampT1))
        algo->setClampT1(1e-6, *clampT1);
    algo->setSequence(spgrSequence);
    auto apply = QI::ApplyF::New();
    apply->SetVerbose(*verbose);
    apply->SetAlgorithm(algo);
    apply->SetOutputAllResiduals(*all_residuals);
    apply->SetPoolsize(*num_threads);
    apply->SetSplitsPerThread(*num_threads); // Unbalanced algorithm
    apply->SetInput(0, data);
    apply->SetMask(*mask);
    apply->SetConst(0, *B1);
    if (*verbose) {
        cout << "Processing" << endl;
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    if (*verbose) {
        cout << "Elapsed time was " << apply->GetTotalTime() << "s" << endl;
        cout << "Writing results files." << endl;
    }
    *outPrefix += "D1_";
    QI::WriteImage(apply->GetOutput(0),*outPrefix + "PD" + QI::OutExt());
    QI::WriteImage(apply->GetOutput(1), *outPrefix + "T1" + QI::OutExt());
    QI::WriteScaledImage(apply->GetResidualOutput(), apply->GetOutput(0), *outPrefix + "residual" + QI::OutExt());
    if (*all_residuals) {
        QI::WriteScaledVectorImage(apply->GetAllResidualsOutput(), apply->GetOutput(0), *outPrefix + "all_residuals" + QI::OutExt());
    }
    if (algo->getIterations() != D1Algo::DefaultIterations) {
        QI::WriteImage(apply->GetIterationsOutput(), *outPrefix + "iterations" + QI::OutExt());
    }
    if (*verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}
