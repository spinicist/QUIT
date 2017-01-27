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
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/NumericalDiff>

#include "cppoptlib/boundedproblem.h"
#include "cppoptlib/solver/lbfgsbsolver.h"
#include "Filters/ImageToVectorFilter.h"
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
    virtual void apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
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
        outputs[0] = QI::clamp(b[1] / (1. - b[0]), m_loPD, m_hiPD);
        outputs[1] = QI::clamp(-m_sequence->TR() / log(b[0]), m_loT1, m_hiT1);
        ArrayXd theory = QI::One_SPGR(m_sequence->flip(), m_sequence->TR(), outputs[0], outputs[1], B1).array().abs();
        ArrayXf r = (data.array() - theory).cast<float>();
        residual = sqrt(r.square().sum() / r.rows());
        resids = itk::VariableLengthVector<float>(r.data(), r.rows());
        its = 1;
    }
};

class D1WLLS : public D1Algo {
public:
    virtual void apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
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
        outputs[0] = QI::clamp(out[0], m_loPD, m_hiPD);
        outputs[1] = QI::clamp(out[1], m_loT1, m_hiT1);
        ArrayXd theory = QI::One_SPGR(m_sequence->flip(), m_sequence->TR(), outputs[0], outputs[1], B1).array().abs();
        ArrayXf r = (data.array() - theory).cast<float>();
        residual = sqrt(r.square().sum() / r.rows());
        resids = itk::VariableLengthVector<float>(r.data(), r.rows());
    }
};

class T1Functor : public DenseFunctor<double> {
    protected:
        const shared_ptr<QI::SequenceBase> m_sequence;
        const ArrayXd m_data;
        const double m_B1;
        const shared_ptr<QI::SCD> m_model;

    public:
        T1Functor(const shared_ptr<QI::SequenceBase> cs, const ArrayXd &data, const double B1) :
            DenseFunctor<double>(2, cs->size()),
            m_sequence(cs), m_data(data), m_B1(B1)
        {
            assert(static_cast<size_t>(m_data.rows()) == values());
        }

        int operator()(const Ref<VectorXd> &p, Ref<ArrayXd> diffs) const {
            eigen_assert(diffs.size() == values());
            ArrayXd s = QI::One_SPGR(m_sequence->flip(), m_sequence->TR(), p[0], p[1], m_B1).array().abs();
            diffs = s - m_data;
            return 0;
        }
};

class D1CostFunction : public cppoptlib::BoundedProblem<double, 2> {
public:
    using BoundedProblem<double, 2>::BoundedProblem;
    using typename Problem<double, 2>::TVector;
    
    shared_ptr<QI::SequenceBase> m_sequence;
    ArrayXd m_data;
    double m_B1;

    Eigen::ArrayXd residuals(const TVector &p) const {
        ArrayXd s = QI::One_SPGR_Magnitude(m_sequence->flip(), m_sequence->TR(), p[0], p[1], m_B1);
        Eigen::ArrayXd diff = s - m_data;
        return diff;
    }

    double value(const TVector &p) {
        return residuals(p).square().sum();
    }
    
    void gradient(const TVector &p, TVector &grad) const {
        ArrayXXd deriv = QI::One_SPGR_Magnitude_Derivs(m_sequence->flip(), m_sequence->TR(), p[0], p[1], m_B1);
        grad = 2*(deriv.colwise()*(residuals(p))).colwise().sum();
    }
};

class D1LBFGSB : public D1Algo {
public:
    virtual void apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                  std::vector<TOutput> &outputs, TConst &residual,
                  TInput &resids, TIters &its) const override
    {
        Eigen::Map<const Eigen::ArrayXf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
        double B1 = consts[0];
        // Improve scaling by dividing the PD down to something sensible.
        // This gets scaled back up at the end.
        const ArrayXd data = indata.cast<double>() / indata.maxCoeff();
        auto stop = cppoptlib::Criteria<double>::defaults();
        stop.iterations = 100;
        stop.gradNorm = 1e-8;
        
        cppoptlib::LbfgsbSolver<D1CostFunction> solver;
        solver.setStopCriteria(stop);
        D1CostFunction cost;
        cost.m_B1 = B1;
        cost.m_data = data;
        cost.m_sequence = this->m_sequence;
        Array2d lower; lower << 0.001, 0.1;
        Array2d upper; upper << 50, 5.0;
        cost.setLowerBound(lower);
        cost.setUpperBound(upper);
        Eigen::Vector2d p; p << 10., 1.0;
        solver.minimize(cost, p);
        outputs[0] = QI::clamp(p[0] * indata.maxCoeff(), m_loPD, m_hiPD);
        outputs[1] = QI::clamp(p[1], m_loT1, m_hiT1);
        if (!isfinite(p[1])) {
            cout << "Not finite" << endl;
            cout << indata.transpose() << endl;
            cout << data.transpose() << endl;
            cout << p.transpose() << endl << B1 << " " << indata.abs().maxCoeff() << endl;
            solver.setDebug(cppoptlib::DebugLevel::High);
            Eigen::Vector2d p; p << 10., 1.0;
            solver.minimize(cost, p);
        }
        ArrayXf r = (cost.residuals(p) * indata.abs().maxCoeff()).cast<float>();
        residual = sqrt(r.square().sum() / r.rows());
        resids = itk::VariableLengthVector<float>(r.data(), r.rows());
    }
};

class D1NLLS : public D1Algo {
public:
    virtual void apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                       std::vector<TOutput> &outputs, TConst &residual,
                       TInput &resids, TIters &its) const override
    {
        Eigen::Map<const Eigen::ArrayXf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
        double B1 = consts[0];
        const ArrayXd data = indata.cast<double>() / indata.maxCoeff();
        T1Functor f(m_sequence, data, B1);
        NumericalDiff<T1Functor> nDiff(f);
        LevenbergMarquardt<NumericalDiff<T1Functor>> lm(nDiff);
        lm.setMaxfev(m_iterations * (m_sequence->size() + 1));
        // PD & T1 - Initial guess of 1s
        VectorXd p(2); p << data.maxCoeff() * 10., 1.;
        lm.minimize(p);
        outputs[0] = QI::clamp(p[0] * indata.maxCoeff(), m_loPD, m_hiPD);
        outputs[1] = QI::clamp(p[1], m_loT1, m_hiT1);
        ArrayXd theory = QI::One_SPGR(m_sequence->flip(), m_sequence->TR(), outputs[0], outputs[1], B1).array().abs();
        ArrayXf r = (data.array() - theory).cast<float>();
        residual = sqrt(r.square().sum() / r.rows());
        resids = itk::VariableLengthVector<float>(r.data(), r.rows());
        its = lm.iterations();
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
    QI::Option<float> clampPD(std::numeric_limits<float>::infinity(),'p',"clampPD","Clamp PD between 0 and value", opts);
    QI::Option<float> clampT1(std::numeric_limits<float>::infinity(),'t',"clampT1","Clamp T1 between 0 and value", opts);
    QI::ImageOption<QI::VolumeF> mask('m', "mask", "Mask input with specified file", opts);
    QI::ImageOption<QI::VolumeF> B1('b', "B1", "B1 Map file (ratio)", opts);
    QI::EnumOption algorithm("lwnb",'l','a',"algo","Choose algorithm (l/w/n/b)", opts);
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
        case 'b': algo = make_shared<D1LBFGSB>(); if (*verbose) cout << "LBFGSB algorithm selected." << endl; break;
    }
    algo->setIterations(*its);
    if (isfinite(*clampPD))
        algo->setClampPD(0, *clampPD);
    if (isfinite(*clampT1))
        algo->setClampT1(0, *clampT1);
    algo->setSequence(spgrSequence);
    auto apply = QI::ApplyF::New();
    apply->SetVerbose(*verbose);
    apply->SetAlgorithm(algo);
    apply->SetOutputAllResiduals(*all_residuals);
    apply->SetPoolsize(*num_threads);
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
    QI::WriteImage(apply->GetOutput(0),*outPrefix + "PD.nii");
    QI::WriteImage(apply->GetOutput(1), *outPrefix + "T1.nii");
    QI::WriteScaledImage(apply->GetResidualOutput(), apply->GetOutput(0), *outPrefix + "residual.nii");
    if (*all_residuals) {
        QI::WriteScaledVectorImage(apply->GetAllResidualsOutput(), apply->GetOutput(0), *outPrefix + "all_residuals.nii");
    }
    if (algo->getIterations() != D1Algo::DefaultIterations) {
        QI::WriteImage(apply->GetIterationsOutput(), *outPrefix + "iterations.nii");
    }
    if (*verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}
