/*
 *  apply_main.cpp
 *
 *  Created by Tobias Wood on 23/01/2012.
 *  Copyright (c) 2012-2013 Tobias Wood.
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

#include "Models.h"
#include "Sequences.h"
#include "Util.h"
#include "Args.h"
#include "IO.h"

//******************************************************************************
// Algorithm Subclasses
//******************************************************************************
class D2Algo : public QI::ApplyF::Algorithm {
protected:
    const std::shared_ptr<QI::SCD> m_model = std::make_shared<QI::SCD>();
    std::shared_ptr<QI::SteadyState> m_sequence;
    size_t m_iterations = 15;
    bool m_elliptical = false;
    double m_loPD = -std::numeric_limits<double>::infinity();
    double m_hiPD = std::numeric_limits<double>::infinity();
    double m_loT2 = -std::numeric_limits<double>::infinity();
    double m_hiT2 = std::numeric_limits<double>::infinity();

public:
    void setIterations(size_t n) { m_iterations = n; }
    void setSequence(std::shared_ptr<QI::SteadyState> &s) { m_sequence = s; }
    void setElliptical(bool e) { m_elliptical = e; }
    void setClampT2(double lo, double hi) { m_loT2 = lo; m_hiT2 = hi; }
    void setClampPD(double lo, double hi) { m_loPD = lo; m_hiPD = hi; }
    size_t numInputs() const override { return m_sequence->count(); }
    size_t numConsts() const override { return 2; }  // T1, B1
    size_t numOutputs() const override { return 2; } // PD, T2
    size_t dataSize() const override { return m_sequence->size(); }
    float zero() const override { return 0.; }

    std::vector<float> defaultConsts() const override {
        std::vector<float> def(2, 1.0f); // T1, B1
        return def;
    }
};

class D2LLS : public D2Algo {
public:
    bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
               std::vector<TOutput> &outputs, TConst &residual,
               TInput &resids, TIters &its) const override
    {
        const double T1 = consts[0];
        const double B1 = consts[1];
        const double TR = m_sequence->TR();
        const double E1 = exp(-TR / T1);
        double PD, T2, E2;
        const Eigen::ArrayXd angles = (m_sequence->flip() * B1);
        Eigen::Map<const Eigen::ArrayXf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
        Eigen::ArrayXd data = indata.cast<double>();
        Eigen::VectorXd Y = data / angles.sin();
        Eigen::MatrixXd X(Y.rows(), 2);
        X.col(0) = data / angles.tan();
        X.col(1).setOnes();
        Eigen::VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
        if (m_elliptical) {
            T2 = 2. * TR / log((b[0]*E1 - 1.) / (b[0] - E1));
            E2 = exp(-TR / T2);
            PD = b[1] * (1. - E1*E2*E2) / (sqrt(E2) * (1. - E1));
        } else {
            T2 = TR / log((b[0]*E1 - 1.)/(b[0] - E1));
            E2 = exp(-TR / T2);
            PD = b[1] * (1. - E1*E2) / (1. - E1);
        }
        Eigen::VectorXd p(5); p << PD, T1, T2, 0, B1;
        Eigen::ArrayXd theory = m_sequence->signal(m_model, p).abs();
        Eigen::ArrayXf r = (data.array() - theory).cast<float>();
        residual = sqrt(r.square().sum() / r.rows());
        resids = itk::VariableLengthVector<float>(r.data(), r.rows());
        outputs[0] = QI::Clamp(PD, m_loPD, m_hiPD);
        outputs[1] = QI::Clamp(T2, m_loT2, m_hiT2);
        its = 1;
        return true;
    }
};

class D2WLLS : public D2Algo {
public:
    bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
               std::vector<TOutput> &outputs, TConst &residual,
               TInput &resids, TIters &its) const override
    {
        const double T1 = consts[0];
        const double B1 = consts[1];
        const double TR = m_sequence->TR();
        const double E1 = exp(-TR / T1);
        double PD, T2, E2;
        const Eigen::ArrayXd angles = (m_sequence->flip() * B1);
        Eigen::Map<const Eigen::ArrayXf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
        Eigen::ArrayXd data = indata.cast<double>();
        Eigen::VectorXd Y = data / angles.sin();
        Eigen::MatrixXd X(Y.rows(), 2);
        X.col(0) = data / angles.tan();
        X.col(1).setOnes();
        Eigen::VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
        if (m_elliptical) {
            T2 = 2. * TR / log((b[0]*E1 - 1.) / (b[0] - E1));
            E2 = exp(-TR / T2);
            PD = b[1] * (1. - E1*E2*E2) / (sqrt(E2) * (1. - E1));
        } else {
            T2 = TR / log((b[0]*E1 - 1.)/(b[0] - E1));
            E2 = exp(-TR / T2);
            PD = b[1] * (1. - E1*E2) / (1. - E1);
        }
        Eigen::VectorXd W(m_sequence->size());
        for (size_t n = 0; n < m_iterations; n++) {
            if (m_elliptical) {
                W = ((1. - E1*E2) * angles.sin() / (1. - E1*E2*E2 - (E1 - E2*E2)*angles.cos())).square();
            } else {
                W = ((1. - E1*E2) * angles.sin() / (1. - E1*E2 - (E1 - E2)*angles.cos())).square();
            }
            b = (X.transpose() * W.asDiagonal() * X).partialPivLu().solve(X.transpose() * W.asDiagonal() * Y);
            if (m_elliptical) {
                T2 = 2. * TR / log((b[0]*E1 - 1.) / (b[0] - E1));
                E2 = exp(-TR / T2);
                PD = b[1] * (1. - E1*E2*E2) / (sqrt(E2) * (1. - E1));
            } else {
                T2 = TR / log((b[0]*E1 - 1.)/(b[0] - E1));
                E2 = exp(-TR / T2);
                PD = b[1] * (1. - E1*E2) / (1. - E1);
            }
        }
        Eigen::VectorXd p(5); p << PD, T1, T2, 0, B1;
        Eigen::ArrayXd theory = m_sequence->signal(m_model, p).abs();
        Eigen::ArrayXf r = (data.array() - theory).cast<float>();
        residual = sqrt(r.square().sum() / r.rows());
        resids = itk::VariableLengthVector<float>(r.data(), r.rows());
        outputs[0] = QI::Clamp(PD, m_loPD, m_hiPD);
        outputs[1] = QI::Clamp(T2, m_loT2, m_hiT2);
        its = m_iterations;
        return true;
    }
};

//******************************************************************************
// T2 Only Functor
//******************************************************************************
class D2Functor : public Eigen::DenseFunctor<double> {
    public:
        const std::shared_ptr<QI::SequenceBase> m_sequence;
        const double m_T1, m_B1;
        const std::shared_ptr<QI::SCD> m_model = std::make_shared<QI::SCD>();
        const Eigen::ArrayXd m_data;

        D2Functor(const double T1, const std::shared_ptr<QI::SequenceBase> s, const Eigen::ArrayXd &d, const double B1, const bool fitComplex, const bool debug = false) :
            DenseFunctor<double>(3, s->size()),
            m_sequence(s), m_data(d),
            m_T1(T1), m_B1(B1)
        {
            assert(static_cast<size_t>(m_data.rows()) == values());
        }

        int operator()(const Eigen::Ref<Eigen::VectorXd> &params, Eigen::Ref<Eigen::ArrayXd> diffs) const {
            eigen_assert(diffs.size() == values());
            Eigen::ArrayXd fullparams(5);
            fullparams << params(0), m_T1, params(1), params(2), m_B1;
            Eigen::ArrayXcd s = m_sequence->signal(m_model, fullparams);
            diffs = s.abs() - m_data;
            return 0;
        }
};

class D2NLLS : public D2Algo {
public:
    bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
               std::vector<TOutput> &outputs, TConst &residual,
               TInput &resids, TIters &its) const override
    {
        const double T1 = consts[0];
        const double B1 = consts[1];
        Eigen::Map<const Eigen::ArrayXf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
        Eigen::ArrayXd data = indata.cast<double>();
        D2Functor f(T1, m_sequence, data, B1, false, false);
        Eigen::NumericalDiff<D2Functor> nDiff(f);
        Eigen::LevenbergMarquardt<Eigen::NumericalDiff<D2Functor>> lm(nDiff);
        lm.setMaxfev(m_iterations * (m_sequence->size() + 1));
        Eigen::VectorXd p(2); p << data.array().maxCoeff() * 5., 0.1;
        lm.minimize(p);
        outputs[0] = QI::Clamp(p[0], m_loPD, m_hiPD);
        outputs[1] = QI::Clamp(p[1], m_loT2, m_hiT2);
        Eigen::VectorXd fullp(5); fullp << outputs[0], T1, outputs[1], 0, B1; // Assume on-resonance
        Eigen::ArrayXd theory = m_sequence->signal(m_model, fullp).abs(); // Sequence will already be elliptical if necessary
        Eigen::ArrayXf r = (data.array() - theory).cast<float>();
        residual = sqrt(r.square().sum() / r.rows());
        resids = itk::VariableLengthVector<float>(r.data(), r.rows());
        its = lm.iterations();
        return true;
    }
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates T2 maps from SSFP data/nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> t1_path(parser, "T1 MAP", "Path to T1 map");
    args::Positional<std::string> ssfp_path(parser, "SSFP FILE", "Path to SSFP data");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag     noprompt(parser, "NOPROMPT", "Suppress input prompts", {'n', "no-prompt"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio) file", {'b', "B1"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> subregion(parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    args::Flag resids(parser, "RESIDS", "Write out residuals for each data-point", {'r', "resids"});
    args::ValueFlag<char> algorithm(parser, "ALGO", "Choose algorithm (l/w/n)", {'a',"algo"}, 'l');
    args::Flag ellipse(parser, "ELLIPTICAL", "Data is band-free ellipse / geometric solution", {'e',"ellipse"});
    args::ValueFlag<int> its(parser, "ITERS", "Max iterations for WLLS/NLLS (default 15)", {'i',"its"}, 15);
    args::ValueFlag<float> clampPD(parser, "CLAMP PD", "Clamp PD between 0 and value", {'p',"clampPD"}, std::numeric_limits<float>::infinity());
    args::ValueFlag<float> clampT2(parser, "CLAMP T2", "Clamp T2 between 0 and value", {'t',"clampT2"}, std::numeric_limits<float>::infinity());
    QI::ParseArgs(parser, argc, argv);
    bool prompt = !noprompt;

    std::shared_ptr<D2Algo> algo;
    switch (algorithm.Get()) {
        case 'l': algo = std::make_shared<D2LLS>();  if (verbose) std::cout << "LLS algorithm selected." << std::endl; break;
        case 'w': algo = std::make_shared<D2WLLS>(); if (verbose) std::cout << "WLLS algorithm selected." << std::endl; break;
        case 'n': algo = std::make_shared<D2NLLS>(); if (verbose) std::cout << "NLLS algorithm selected." << std::endl; break;
    }

    algo->setIterations(its);
    if (clampPD) algo->setClampPD(0, clampPD.Get());
    if (clampT2) algo->setClampT2(0, clampT2.Get());
    std::shared_ptr<QI::SteadyState> ssfp;
    if (ellipse) {
        ssfp = std::make_shared<QI::SSFP_GS>(std::cin, prompt);
    } else {
        ssfp = std::make_shared<QI::SSFPSimple>(std::cin, prompt);
    }
    if (verbose) std::cout << *ssfp << std::endl;
    algo->setSequence(ssfp);
    algo->setElliptical(ellipse);

    if (verbose) std::cout << "Reading T1 Map from: " << QI::CheckPos(t1_path) << std::endl;
    auto T1 = QI::ReadImage(QI::CheckPos(t1_path));

    if (verbose) std::cout << "Opening SSFP file: " << QI::CheckPos(ssfp_path) << std::endl;
    auto data = QI::ReadVectorImage<float>(QI::CheckPos(ssfp_path));
    auto apply = QI::ApplyF::New();
    apply->SetAlgorithm(algo);
    apply->SetOutputAllResiduals(resids);
    apply->SetPoolsize(threads.Get());
    apply->SetInput(0, data);
    apply->SetConst(0, T1);
    if (B1) apply->SetConst(1, QI::ReadImage(B1.Get()));
    if (mask) apply->SetMask(QI::ReadImage(mask.Get()));

    if (verbose) {
        std::cout << "apply setup complete. Processing." << std::endl;
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    if (verbose) {
        std::cout << "Elapsed time was " << apply->GetTotalTime() << "s" << std::endl;
        std::cout << "Writing results files." << std::endl;

    }
    std::string outPrefix = outarg.Get() + "D2_";
    QI::WriteImage(apply->GetOutput(0), outPrefix + "PD" + QI::OutExt());
    QI::WriteImage(apply->GetOutput(1), outPrefix + "T2" + QI::OutExt());
    QI::WriteScaledImage(apply->GetResidualOutput(), apply->GetOutput(0), outPrefix + "residual" + QI::OutExt());
    if (resids) {
        QI::WriteScaledVectorImage(apply->GetAllResidualsOutput(), apply->GetOutput(0), outPrefix + "all_residuals" + QI::OutExt());
    }
    if (verbose) std::cout << "All done." << std::endl;
    return EXIT_SUCCESS;
}
