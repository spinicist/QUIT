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

#include "cppoptlib/meta.h"
#include "cppoptlib/boundedproblem.h"
#include "cppoptlib/solver/lbfgsbsolver.h"
#include "cppoptlib/solver/cmaesbsolver.h"
#include "QI/Util.h"
#include "QI/Option.h"
#include "QI/Models/Models.h"
#include "QI/Sequences/Sequences.h"
#include "Filters/ApplyAlgorithmFilter.h"

using namespace std;
using namespace Eigen;

class FMAlgo : public QI::ApplyF::Algorithm {
protected:
	shared_ptr<QI::SSFPSimple> m_sequence;
    bool m_asymmetric = false;

public:
    void setSequence(shared_ptr<QI::SSFPSimple> s) { m_sequence = s; }
    void setAsymmetric(const bool b) { m_asymmetric = b; }
    
    size_t numInputs() const override  { return m_sequence->count(); }
    size_t numConsts() const override  { return 2; }
    size_t numOutputs() const override { return 3; }
    size_t dataSize() const override   { return m_sequence->size(); }
    const float &zero(const size_t i) const override { static float zero = 0; return zero; }
    virtual std::vector<float> defaultConsts() const override {
        std::vector<float> def(2, 1.0f); // T1 & B1
        return def;
    }
};

class FMCostFunction : public cppoptlib::BoundedProblem<double, 3> {
public:
    using BoundedProblem<double, 3>::BoundedProblem;
    using typename Problem<double, 3>::TVector;
    
    Eigen::ArrayXd m_data;
    double m_T1, m_B1;
    shared_ptr<QI::SSFPSimple> m_sequence;
    const shared_ptr<QI::SCD> m_model = make_shared<QI::SCD>();

    Eigen::ArrayXd residuals(const Eigen::VectorXd &p) const {
        ArrayXd s = QI::One_SSFP_Echo_Magnitude(m_sequence->flip(), m_sequence->phase_incs(), m_sequence->TR(), p[0], m_T1, p[1], p[2]/m_sequence->TR(), m_B1);
        Eigen::ArrayXd diff = s - m_data;
        return diff;
    }

    double value(const TVector &p) {
        return residuals(p).square().sum();
    }
    
    void gradient(const TVector &p, TVector &grad) const {
        ArrayXXd deriv = QI::One_SSFP_Echo_Derivs(m_sequence->flip(), m_sequence->phase_incs(), m_sequence->TR(), p[0], m_T1, p[1], p[2]/m_sequence->TR(), m_B1);
        grad = 2*(deriv.colwise()*(residuals(p))).colwise().sum();
    }
};

class BFGSAlgo : public FMAlgo {
protected:
    int m_nstart = 2;
public:
    BFGSAlgo(const int starts) : FMAlgo(), m_nstart(starts) {}

    virtual void apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                  std::vector<TOutput> &outputs, TConst &residual,
                  TInput &resids, TIters &its) const override
    {
        const double T1 = consts[0];
        const double B1 = consts[1];
        if (isfinite(T1) && (T1 > 0.001)) {
            // Improve scaling by dividing the PD down to something sensible.
            // This gets scaled back up at the end.
            Eigen::Map<const Eigen::ArrayXf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
            ArrayXd data = indata.cast<double>() / indata.maxCoeff();
            auto stop = cppoptlib::Criteria<double>::defaults();
            stop.iterations = 100;
            stop.gradNorm = 1e-8;
            
            cppoptlib::LbfgsbSolver<FMCostFunction> solver;
            solver.setStopCriteria(stop);
            FMCostFunction cost;
            cost.m_B1 = B1;
            cost.m_data = data;
            cost.m_sequence = this->m_sequence;
            cost.m_T1 = T1;
            Array3d lower; lower << 0., 2.*this->m_sequence->TR(), -0.55;
            Array3d upper; upper << 20,    T1,                         0.55;
            if (!this->m_asymmetric)
                lower[2] = 0.;
            cost.setLowerBound(lower);
            cost.setUpperBound(upper);
            
            // f0 is scaled by TR in the cost function so that scaling is better here
            vector<double> f0_starts(m_nstart);
            if (this->m_asymmetric) {
                double step = (upper[2] - lower[2]) / (m_nstart+1);
                for (int i = 0; i < m_nstart; i++) {
                    f0_starts[i] = lower[2] + (i+1)*step;
                }
            } else {
                for (int i = 0; i < m_nstart; i++) {
                    f0_starts[i] = 0.001 + upper[2]*(static_cast<double>(i) / static_cast<double>(m_nstart));
                }
            }

            double best = numeric_limits<double>::infinity();
            Eigen::Array3d bestP;
            its = 0;
            for (const double &f0 : f0_starts) {
                Eigen::Vector3d p; p << 5., 0.1 * T1, f0; // Yarnykh gives T2 = 0.045 * T1 in brain, but best to overestimate for CSF
                solver.minimize(cost, p);
                double r = cost(p);
                if (r < best) {
                    best = r;
                    bestP = p;
                    its = solver.criteria().iterations;
                }
            }
            outputs[0] = bestP[0] * indata.abs().maxCoeff();
            outputs[1] = bestP[1];
            outputs[2] = bestP[2] / this->m_sequence->TR();
            ArrayXf r = cost.residuals(bestP).cast<float>();
            residual = sqrt(r.square().sum() / r.rows());
            resids = itk::VariableLengthVector<float>(r.data(), r.rows());
        } else {
            outputs[0] = 0.;
            outputs[1] = 0.;
            outputs[2] = 0.;
            residual = 0;
            resids.Fill(0.);
            its = 0;
        }
    }
};

class CMAlgo : public FMAlgo {
public:
    virtual void apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                  std::vector<TOutput> &outputs, TConst &residual,
                  TInput &resids, TIters &its) const override
    {
        const double T1 = consts[0];
        const double B1 = consts[1];
        if (isfinite(T1) && (T1 > 0.001)) {
            // Improve scaling by dividing the PD down to something sensible.
            // This gets scaled back up at the end.
            Eigen::Map<const Eigen::ArrayXf> indata(inputs[0].GetDataPointer(), inputs[0].Size());
            ArrayXd data = indata.cast<double>() / indata.abs().maxCoeff();
            cppoptlib::CMAesBSolver<FMCostFunction> solver;
            // f0 is scaled by TR in the cost function so that scaling is better here
            Array3d lower; lower << 0., 2.*this->m_sequence->TR(), -0.55;
            Array3d upper; upper << 20,    T1,                         0.55;
            if (!this->m_asymmetric)
                lower[2] = 0.;
            FMCostFunction cost(lower, upper);
            cost.m_B1 = B1;
            cost.m_data = data;
            cost.m_sequence = this->m_sequence;
            cost.m_T1 = T1;
            
            Eigen::Vector3d p; p << 5., 0.1 * T1, 0.1; // Yarnykh gives T2 = 0.045 * T1 in brain, but best to overestimate for CSF
            //solver.setDebug(cppoptlib::DebugLevel::Low);
            solver.minimize(cost, p);
            its = solver.criteria().iterations;
            outputs[0] = p[0] * indata.abs().maxCoeff();
            outputs[1] = p[1];
            outputs[2] = p[2] / this->m_sequence->TR();
            ArrayXf rf = cost.residuals(p).cast<float>() * indata.abs().maxCoeff();
            residual = sqrt(rf.square().sum() / rf.rows());
            resids = itk::VariableLengthVector<float>(rf.data(), rf.rows());
        } else {
            outputs[0] = 0.;
            outputs[1] = 0.;
            outputs[2] = 0.;
            residual = 0;
            resids.Fill(0.);
            its = 0;
        }
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
    QI::RegionOption<QI::ApplyF::TRegion> subregion('s',"subregion","Process subregion starting at voxel I,J,K with size SI,SJ,SK", opts);
    QI::Switch flex('f',"flex", "Specify all phase-incs for all flip-angles", opts);
    QI::Option<int> nstart(2, 'F',"off","Number of off-resonance start points", opts);
    QI::Switch asym('A',"asym","Fit +/- off-resonance frequency", opts);
    QI::ImageOption<QI::VolumeF> mask('m', "mask", "Mask input with specified file", opts);
    QI::ImageOption<QI::VolumeF> B1('b', "B1", "B1 Map file (ratio)", opts);
    QI::EnumOption algorithm("bc",'b','a',"algo","Choose algorithm (b/c)", opts);
    QI::Option<std::string> outPrefix("", 'o', "out","Add a prefix to output filenames", opts);
    QI::Switch suppress('n',"no-prompt","Suppress input prompts", opts);
    QI::Switch verbose('v',"verbose","Print more information", opts);
    QI::Help help(opts);
    std::vector<std::string> nonopts = opts.parse(argc, argv);
    if (nonopts.size() != 2) {
        std::cerr << opts << std::endl;
        std::cerr << "Wrong number of arguments. Need a T1 map and one SSFP file." << std::endl;
        return EXIT_FAILURE;
    }

    shared_ptr<QI::SSFPSimple> ssfpSequence;
    if (*flex)
        ssfpSequence = make_shared<QI::SSFPEchoFlex>(cin, !(*suppress));
    else
        ssfpSequence = make_shared<QI::SSFPEcho>(cin, !(*suppress));
    if (*verbose) cout << *ssfpSequence << endl;

    if (*verbose) cout << "Reading T1 Map from: " << nonopts[0] << endl;
    auto T1 = QI::ReadImage(nonopts[0]);
    if (*verbose) cout << "Opening SSFP file: " << nonopts[1] << endl;
    auto ssfpData = QI::ReadVectorImage<float>(nonopts[1]);
    auto apply = QI::ApplyF::New();
    shared_ptr<FMAlgo> algo;
    switch (*algorithm) {
        case 'b': algo = make_shared<BFGSAlgo>(*nstart); if (*verbose) cout << "LBFGSB algorithm selected." << endl; break;
        case 'c': algo = make_shared<CMAlgo>(); if (*verbose) cout << "CM algorithm selected." << endl; break;
        default: throw(std::runtime_error("Invalid algorithm: " + std::to_string(*algorithm)));
    }

    algo->setSequence(ssfpSequence);
    algo->setAsymmetric(*asym);
    apply->SetVerbose(*verbose);
    apply->SetAlgorithm(algo);
    apply->SetOutputAllResiduals(*all_residuals);
    apply->SetPoolsize(*num_threads);
    apply->SetInput(0, ssfpData);
    apply->SetConst(0, T1);
    apply->SetConst(1, *B1);
    apply->SetMask(*mask);
    if (subregion.set()) {
        if (*verbose) cout << "Setting subregion: " << *subregion << endl;
        apply->SetSubregion(*subregion);
    }
    if (*verbose) {
        cout << "Processing" << endl;
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    if (*verbose) {
        cout << "Elapsed time was " << apply->GetTotalTime() << "s" << endl;
        cout << "Mean time per voxel was " << apply->GetMeanTime() << "s" << endl;
        cout << "Writing results files." << endl;
    }
    *outPrefix = *outPrefix + "FM_";
    QI::WriteImage(apply->GetOutput(0), *outPrefix + "PD.nii");
    QI::WriteImage(apply->GetOutput(1), *outPrefix + "T2.nii");
    QI::WriteImage(apply->GetOutput(2), *outPrefix + "f0.nii");
    QI::WriteImage(apply->GetIterationsOutput(), *outPrefix + "its.nii");
    QI::WriteScaledImage(apply->GetResidualOutput(), apply->GetOutput(0), *outPrefix + "residual.nii");
    if (*all_residuals) {
        QI::WriteScaledVectorImage(apply->GetAllResidualsOutput(), apply->GetOutput(0), *outPrefix + "all_residuals.nii");
    }
    return EXIT_SUCCESS;
}
