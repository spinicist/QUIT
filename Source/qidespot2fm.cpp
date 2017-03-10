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
#include "QI/Util.h"
#include "QI/Args.h"
#include "QI/IO.h"
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

class FMCost : public ceres::CostFunction {
private:
    const Eigen::ArrayXd &m_data;
    double m_T1, m_B1;
    shared_ptr<QI::SSFPSimple> m_sequence;

public:
    FMCost(const Eigen::ArrayXd &d, shared_ptr<QI::SSFPSimple> s,
           const double T1, const double B1) :
        m_data(d), m_sequence(s), m_T1(T1), m_B1(B1)
    {
        mutable_parameter_block_sizes()->push_back(3);
        set_num_residuals(d.size());
    }

    Eigen::ArrayXd residuals(const Eigen::VectorXd &p) const {
        ArrayXd s = QI::One_SSFP_Echo_Magnitude(m_sequence->allFlip(), m_sequence->allPhi(), m_sequence->TR(), p[0], m_T1, p[1], p[2], m_B1);
        Eigen::ArrayXd diff = s - m_data;
        return diff;
    }

    virtual bool Evaluate(double const* const* parameters,
                        double* resids,
                        double** jacobians) const override
    {
        Eigen::Map<const Eigen::Array3d> p(parameters[0]);
        Eigen::Map<Eigen::ArrayXd> r(resids, m_data.size());
        r = residuals(p);
        if (jacobians && jacobians[0]) {
            Eigen::Map<Eigen::Matrix<double, -1, -1, RowMajor>> j(jacobians[0], m_data.size(), 3);
            j = QI::One_SSFP_Echo_Derivs(m_sequence->allFlip(), m_sequence->allPhi(), m_sequence->TR(), p[0], m_T1, p[1], p[2], m_B1);
        }
        return true;
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

            vector<double> f0_starts = {0, 0.4/m_sequence->TR()};
            if (this->m_asymmetric)
                f0_starts.push_back(-0.4/m_sequence->TR());

            double best = numeric_limits<double>::infinity();
            Eigen::Array3d p, bestP;
            ceres::Problem problem;
            problem.AddResidualBlock(new FMCost(data, m_sequence, T1, B1), NULL, p.data());
            problem.SetParameterLowerBound(p.data(), 0, 1.);
            problem.SetParameterLowerBound(p.data(), 1, m_sequence->TR());
            problem.SetParameterUpperBound(p.data(), 1, T1);
            problem.SetParameterLowerBound(p.data(), 2, -1./(2.*m_sequence->TR()));
            problem.SetParameterUpperBound(p.data(), 2,  1./(2.*m_sequence->TR()));
            ceres::Solver::Options options;
            ceres::Solver::Summary summary;
            options.max_num_iterations = 50;
            options.function_tolerance = 1e-5;
            options.gradient_tolerance = 1e-6;
            options.parameter_tolerance = 1e-4;
            options.logging_type = ceres::SILENT;
            for (const double &f0 : f0_starts) {
                p = {10., 0.1 * T1, f0}; // Yarnykh gives T2 = 0.045 * T1 in brain, but best to overestimate for CSF
                ceres::Solve(options, &problem, &summary);
                double r = summary.final_cost;
                if (r < best) {
                    best = r;
                    bestP = p;
                    its = summary.iterations.size();
                }
            }
            outputs[0] = bestP[0] * indata.maxCoeff();
            outputs[1] = bestP[1];
            outputs[2] = bestP[2];
            
            residual = best;
            if (resids.Size() > 0) {
                assert(resids.Size() == data.size());
                vector<double> r_temp(data.size());
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
    }
};

/*class CMAlgo : public FMAlgo {
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
};*/

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();

    QI::ArgParser args{argc, argv,
        "Usage is: qidespot2fm t1_map ssfp_input [options]\n"
        "Calculates a T2 map from SSFP data and a T1 map.",
        {{"help",      'h', "Display the help message and quit", false},
         {"verbose",   'v', "Print more information", false},
         {"no-prompt", 'n', "Suppress input prompts", false},
         {"threads",   'T', "Use N threads (default=4, 0=hardware limit)", true},
         {"out",       'o', "Add a prefix to output filenames", true},
         {"rename",    'r', "Rename using specified header field", true},
         {"algo",      'a', "Choose algorithm (b/c)", true},
         {"B1",        'b', "B1 Map file (ratio)", true},
         {"mask",      'm', "Only process within mask file", true},
         {"asym",      'A', "Fit +/- off-resonance frequency", false},
         {"off",       'F', "Number of off-resonance start points", true},
         {"flex",      'f', "Flexible input (do not expand incs)", false},
         {"subregion", 's', "Process subregion starting at voxel I,J,K with size SI,SJ,SK", false},
         {"resids",    'r', "Write out residuals for each data-point", false}}
    };

    bool verbose = args.option_present("verbose");
    bool prompt  = !args.option_present("suppress");

    std::deque<std::string> nonopts = args.nonoptions();
    if (nonopts.size() != 2) {
        std::cerr << "Wrong number of arguments. Need a T1 map and one SSFP file." << std::endl;
        return EXIT_FAILURE;
    }

    shared_ptr<QI::SSFPSimple> ssfpSequence;
    if (args.option_present("flex"))
        ssfpSequence = make_shared<QI::SSFPEchoFlex>(cin, prompt);
    else
        ssfpSequence = make_shared<QI::SSFPEcho>(cin, prompt);
    if (verbose) cout << *ssfpSequence << endl;

    if (verbose) cout << "Reading T1 Map from: " << nonopts[0] << endl;
    auto T1 = QI::ReadImage(nonopts[0]);
    if (verbose) cout << "Opening SSFP file: " << nonopts[1] << endl;
    auto ssfpData = QI::ReadVectorImage<float>(nonopts[1]);
    auto apply = QI::ApplyF::New();
    shared_ptr<FMAlgo> algo;
    switch (args.option_value("algo", 'b')) {
        case 'b': algo = make_shared<BFGSAlgo>(args.option_value("off", 2)); if (verbose) cout << "LBFGSB algorithm selected." << endl; break;
        // case 'c': algo = make_shared<CMAlgo>(); if (verbose) cout << "CM algorithm selected." << endl; break;
        default: throw(std::runtime_error("Invalid algorithm specified"));
    }

    algo->setSequence(ssfpSequence);
    algo->setAsymmetric(args.option_present("asym"));
    apply->SetVerbose(verbose);
    apply->SetAlgorithm(algo);
    apply->SetOutputAllResiduals(args.option_present("resids"));
    const int nthreads = args.option_value("threads", 4);
    std::cout << "Using " << nthreads << " threads" << std::endl;
    apply->SetPoolsize(nthreads);
    apply->SetSplitsPerThread(nthreads); // Fairly unbalanced algorithm
    apply->SetInput(0, ssfpData);
    apply->SetConst(0, T1);
    apply->SetConst(1, QI::ReadOptImage("B1", args));
    apply->SetMask(QI::ReadOptImage("mask", args));
    if (args.option_present("subregion")) {
        apply->SetSubregion(QI::RegionOpt("subregion", args));
    }
    if (verbose) {
        cout << "Processing" << endl;
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    if (verbose) {
        cout << "Elapsed time was " << apply->GetTotalTime() << "s" << endl;
        cout << "Writing results files." << endl;
    }
    std::string outPrefix = args.string_value("out", "") + "FM_";
    QI::WriteImage(apply->GetOutput(0), outPrefix + "PD.nii");
    QI::WriteImage(apply->GetOutput(1), outPrefix + "T2.nii");
    QI::WriteImage(apply->GetOutput(2), outPrefix + "f0.nii");
    QI::WriteImage(apply->GetIterationsOutput(), outPrefix + "its.nii");
    QI::WriteScaledImage(apply->GetResidualOutput(), apply->GetOutput(0), outPrefix + "residual.nii");
    if (args.option_present("resids")) {
        QI::WriteScaledVectorImage(apply->GetAllResidualsOutput(), apply->GetOutput(0), outPrefix + "all_residuals.nii");
    }
    return EXIT_SUCCESS;
}
