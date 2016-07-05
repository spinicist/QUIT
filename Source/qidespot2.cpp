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
 
#include "Filters/ImageToVectorFilter.h"
#include "QI/Models/Models.h"
#include "QI/Sequences/Sequences.h"
#include "QI/Util.h"
#include "QI/Option.h"

using namespace std;
using namespace Eigen;

//******************************************************************************
// Algorithm Subclasses
//******************************************************************************
class D2Algo : public QI::ApplyF::Algorithm {
protected:
    const shared_ptr<QI::SCD> m_model = make_shared<QI::SCD>();
    shared_ptr<QI::SteadyState> m_sequence;
    size_t m_iterations = 15;
    bool m_elliptical = false;
    double m_loPD = -numeric_limits<double>::infinity();
    double m_hiPD = numeric_limits<double>::infinity();
    double m_loT2 = -numeric_limits<double>::infinity();
    double m_hiT2 = numeric_limits<double>::infinity();

public:
    void setIterations(size_t n) { m_iterations = n; }
    void setSequence(shared_ptr<QI::SteadyState> &s) { m_sequence = s; }
    void setElliptical(bool e) { m_elliptical = e; }
    void setClampT2(double lo, double hi) { m_loT2 = lo; m_hiT2 = hi; }
    void setClampPD(double lo, double hi) { m_loPD = lo; m_hiPD = hi; }
    size_t numInputs() const override { return m_sequence->count(); }
    size_t numConsts() const override { return 2; }  // T1, B1
    size_t numOutputs() const override { return 2; } // PD, T2
    size_t dataSize() const override { return m_sequence->size(); }
    const float &zero(const size_t i) const override { static float zero = 0; return zero; }

    virtual std::vector<float> defaultConsts() override {
        std::vector<float> def(2, 1.0f); // T1, B1
        return def;
    }
};

class D2LLS : public D2Algo {
public:
    virtual void apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                  std::vector<TOutput> &outputs, TConst &residual,
                  TInput &resids, TIters &its) const override
    {
        const double T1 = consts[0];
        const double B1 = consts[1];
        const double TR = m_sequence->TR();
        const double E1 = exp(-TR / T1);
        double PD, T2, E2;
        const ArrayXd angles = (m_sequence->flip() * B1);
        Eigen::Map<const Eigen::ArrayXf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
        ArrayXd data = indata.cast<double>();
        VectorXd Y = data / angles.sin();
        MatrixXd X(Y.rows(), 2);
        X.col(0) = data / angles.tan();
        X.col(1).setOnes();
        VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
        if (m_elliptical) {
            T2 = 2. * TR / log((b[0]*E1 - 1.) / (b[0] - E1));
            E2 = exp(-TR / T2);
            PD = b[1] * (1. - E1*E2*E2) / (sqrt(E2) * (1. - E1));
        } else {
            T2 = TR / log((b[0]*E1 - 1.)/(b[0] - E1));
            E2 = exp(-TR / T2);
            PD = b[1] * (1. - E1*E2) / (1. - E1);
        }
        VectorXd p(5); p << PD, T1, T2, 0, B1;
        ArrayXd theory = m_sequence->signal(m_model, p).abs();
        ArrayXf r = (data.array() - theory).cast<float>();
        residual = sqrt(r.square().sum() / r.rows());
        resids = itk::VariableLengthVector<float>(r.data(), r.rows());
        outputs[0] = QI::clamp(PD, m_loPD, m_hiPD);
        outputs[1] = QI::clamp(T2, m_loT2, m_hiT2);
        its = 1;
    }
};

class D2WLLS : public D2Algo {
public:
    virtual void apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                  std::vector<TOutput> &outputs, TConst &residual,
                  TInput &resids, TIters &its) const override
    {
        const double T1 = consts[0];
        const double B1 = consts[1];
        const double TR = m_sequence->TR();
        const double E1 = exp(-TR / T1);
        double PD, T2, E2;
        const ArrayXd angles = (m_sequence->flip() * B1);
        Eigen::Map<const Eigen::ArrayXf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
        ArrayXd data = indata.cast<double>();
        VectorXd Y = data / angles.sin();
        MatrixXd X(Y.rows(), 2);
        X.col(0) = data / angles.tan();
        X.col(1).setOnes();
        VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
        if (m_elliptical) {
            T2 = 2. * TR / log((b[0]*E1 - 1.) / (b[0] - E1));
            E2 = exp(-TR / T2);
            PD = b[1] * (1. - E1*E2*E2) / (sqrt(E2) * (1. - E1));
        } else {
            T2 = TR / log((b[0]*E1 - 1.)/(b[0] - E1));
            E2 = exp(-TR / T2);
            PD = b[1] * (1. - E1*E2) / (1. - E1);
        }
        VectorXd W(m_sequence->size());
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
        VectorXd p(5); p << PD, T1, T2, 0, B1;
        ArrayXd theory = m_sequence->signal(m_model, p).abs();
        ArrayXf r = (data.array() - theory).cast<float>();
        residual = sqrt(r.square().sum() / r.rows());
        resids = itk::VariableLengthVector<float>(r.data(), r.rows());
        outputs[0] = QI::clamp(PD, m_loPD, m_hiPD);
        outputs[1] = QI::clamp(T2, m_loT2, m_hiT2);
        its = m_iterations;
    }
};

//******************************************************************************
// T2 Only Functor
//******************************************************************************
class D2Functor : public DenseFunctor<double> {
    public:
        const shared_ptr<QI::SequenceBase> m_sequence;
        const double m_T1, m_B1;
        const shared_ptr<QI::SCD> m_model = make_shared<QI::SCD>();
        const ArrayXd m_data;

        D2Functor(const double T1, const shared_ptr<QI::SequenceBase> s, const ArrayXd &d, const double B1, const bool fitComplex, const bool debug = false) :
            DenseFunctor<double>(3, s->size()),
            m_sequence(s), m_data(d),
            m_T1(T1), m_B1(B1)
        {
            assert(static_cast<size_t>(m_data.rows()) == values());
        }

        int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
            eigen_assert(diffs.size() == values());
            ArrayXd fullparams(5);
            fullparams << params(0), m_T1, params(1), params(2), m_B1;
            ArrayXcd s = m_sequence->signal(m_model, fullparams);
            diffs = s.abs() - m_data;
            return 0;
        }
};

class D2NLLS : public D2Algo {
public:
    virtual void apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                  std::vector<TOutput> &outputs, TConst &residual,
                  TInput &resids, TIters &its) const override
    {
        const double T1 = consts[0];
        const double B1 = consts[1];
        Eigen::Map<const Eigen::ArrayXf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
        ArrayXd data = indata.cast<double>();
        D2Functor f(T1, m_sequence, data, B1, false, false);
        NumericalDiff<D2Functor> nDiff(f);
        LevenbergMarquardt<NumericalDiff<D2Functor>> lm(nDiff);
        lm.setMaxfev(m_iterations * (m_sequence->size() + 1));
        VectorXd p(2); p << data.array().maxCoeff() * 5., 0.1;
        lm.minimize(p);
        outputs[0] = QI::clamp(p[0], m_loPD, m_hiPD);
        outputs[1] = QI::clamp(p[1], m_loT2, m_hiT2);
        VectorXd fullp(5); fullp << outputs[0], T1, outputs[1], 0, B1; // Assume on-resonance
        ArrayXd theory = m_sequence->signal(m_model, fullp).abs(); // Sequence will already be elliptical if necessary
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
    QI::OptionList opts("Usage is: qidespot2 [options] T1_map SSFP_file");
    QI::Switch all_residuals('r',"resids","Write out per flip-angle residuals", opts);
    QI::Option<int> num_threads(4,'T',"threads","Use N threads (default=4, 0=hardware limit)", opts);
    QI::Option<int> its(15,'i',"its","Max iterations for WLLS/NLLS (default 15)", opts);
    QI::Option<float> clampPD(std::numeric_limits<float>::infinity(),'p',"clampPD","Clamp PD between 0 and value", opts);
    QI::Option<float> clampT2(std::numeric_limits<float>::infinity(),'t',"clampT2","Clamp T2 between 0 and value", opts);
    QI::Option<std::string> outPrefix("", 'o', "out","Add a outPrefix to output filenames", opts);
    QI::ImageOption<QI::VolumeF> mask('m', "mask", "Mask input with specified file", opts);
    QI::ImageOption<QI::VolumeF> B1('b', "B1", "B1 Map file (ratio)", opts);
    QI::EnumOption algorithm("lwnb",'l','a',"algo","Choose algorithm (l/w/n)", opts);
    QI::Switch elliptical('e',"elliptical","Input is band-free elliptical/GS data", opts);
    QI::Switch suppress('n',"no-prompt","Suppress input prompts", opts);
    QI::Switch verbose('v',"verbose","Print more information", opts);
    QI::Help help(opts);
    std::vector<std::string> nonopts = opts.parse(argc, argv);
    if (nonopts.size() != 2) {
        std::cerr << opts << std::endl;
        std::cerr << "T1 map and SSFP filename required as inputs." << std::endl;
        return EXIT_FAILURE;
    }

    shared_ptr<D2Algo> algo;
    switch (*algorithm) {
        case 'l': algo = make_shared<D2LLS>();  if (*verbose) cout << "LLS algorithm selected." << endl; break;
        case 'w': algo = make_shared<D2WLLS>(); if (*verbose) cout << "WLLS algorithm selected." << endl; break;
        case 'n': algo = make_shared<D2NLLS>(); if (*verbose) cout << "NLLS algorithm selected." << endl; break;
    }

    algo->setIterations(*its);
    if (isfinite(*clampPD))
        algo->setClampPD(0, *clampPD);
    if (isfinite(*clampT2))
        algo->setClampT2(0, *clampT2);
    shared_ptr<QI::SteadyState> ssfp;
    if (*elliptical) {
        ssfp = make_shared<QI::SSFP_GS>(cin, !*suppress);
    } else {
        ssfp = make_shared<QI::SSFPSimple>(cin, !*suppress);
    }
    if (*verbose) cout << *ssfp << endl;
    algo->setSequence(ssfp);
    algo->setElliptical(*elliptical);

    if (*verbose) cout << "Reading T1 Map from: " << nonopts[0] << endl;
    auto T1 = QI::ReadImage(nonopts[0]);

    if (*verbose) cout << "Opening SSFP file: " << nonopts[1] << endl;
    auto data = QI::ReadVectorImage<float>(nonopts[1]);
    auto apply = itk::ApplyAlgorithmFilter<D2Algo>::New();
    apply->SetAlgorithm(algo);
    apply->SetOutputAllResiduals(*all_residuals);
    apply->SetPoolsize(*num_threads);
    apply->SetInput(0, data);
    apply->SetConst(0, T1);
    apply->SetConst(1, *B1);
    apply->SetMask(*mask);

    if (*verbose) {
        cout << "apply setup complete. Processing." << endl;
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    if (*verbose) {
        cout << "Elapsed time was " << apply->GetTotalTime() << "s" << endl;
        cout << "Mean time per voxel was " << apply->GetMeanTime() << "s" << endl;
        cout << "Writing results files." << endl;

    }
    *outPrefix = *outPrefix + "D2_";
    QI::WriteImage(apply->GetOutput(0), *outPrefix + "PD.nii");
    QI::WriteImage(apply->GetOutput(1), *outPrefix + "T2.nii");
    QI::WriteScaledImage(apply->GetResidualOutput(), apply->GetOutput(0), *outPrefix + "residual.nii");
    if (*all_residuals) {
        QI::WriteScaledVectorImage(apply->GetAllResidualsOutput(), apply->GetOutput(0), *outPrefix + "all_residuals.nii");
    }
    if (*verbose) cout << "All done." << endl;
    return EXIT_SUCCESS;
}
