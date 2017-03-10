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
#include "args.hxx"
#include "ceres/ceres.h"

#include "QI/Util.h"
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
            Eigen::Map<Eigen::Matrix<double, -1, -1, RowMajor>> j(jacobians[0], m_data.size(), p.size());
            j = QI::One_SSFP_Echo_Derivs(m_sequence->allFlip(), m_sequence->allPhi(), m_sequence->TR(), p[0], m_T1, p[1], p[2], m_B1);
        }
        return true;
    }

};

class LM_FM : public FMAlgo {
protected:

public:
    LM_FM() : FMAlgo() {}

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
            
            residual = best * indata.maxCoeff();
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

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Calculates a T2 map from SSFP data and a T1 map. Usage is: qidespot2fm t1_map ssfp_input [options]", "tobias.wood@kcl.ac.uk");
    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag     prompt(parser, "PROMPT", "Suppress input prompts", {'n', "no-prompt"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> B1(parser, "B1", "B1 map (ratio) file", {'b', "B1"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::Flag asym(parser, "ASYM", "Fit +/- off-resonance frequency", {'A', "asym"});
    args::Flag flex(parser, "FLEX", "Flexible input (do not tile flip-angles/phase-incs)", {'f', "flex"});
    args::ValueFlag<std::string> subregion(parser, "SUBREGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    args::Flag resids(parser, "RESIDS", "Write out residuals for each data-point", {'r', "resids"});
    args::PositionalList<std::string> files(parser, "FILES", "Input files - T1 map and SSFP");
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return EXIT_SUCCESS;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return EXIT_FAILURE;
    } catch (args::ValidationError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return EXIT_FAILURE;
    }
    auto filelist = args::get(files);
    if (filelist.size() != 2) {
        std::cerr << "Wrong number of positional arguments." << std::endl;
        return EXIT_FAILURE;
    }

    shared_ptr<QI::SSFPSimple> ssfpSequence;
    if (flex)
        ssfpSequence = make_shared<QI::SSFPEchoFlex>(cin, prompt);
    else
        ssfpSequence = make_shared<QI::SSFPEcho>(cin, prompt);
    if (verbose) cout << *ssfpSequence << endl;

    if (verbose) cout << "Reading T1 Map from: " << filelist[0] << endl;
    auto T1 = QI::ReadImage(filelist[0]);
    if (verbose) cout << "Opening SSFP file: " << filelist[1] << endl;
    auto ssfpData = QI::ReadVectorImage<float>(filelist[1]);
    auto apply = QI::ApplyF::New();
    shared_ptr<FMAlgo> algo = make_shared<LM_FM>();

    algo->setSequence(ssfpSequence);
    algo->setAsymmetric(asym);
    apply->SetVerbose(verbose);
    apply->SetAlgorithm(algo);
    apply->SetOutputAllResiduals(resids);
    if (verbose) std::cout << "Using " << args::get(threads) << " threads" << std::endl;
    apply->SetPoolsize(args::get(threads));
    apply->SetSplitsPerThread(args::get(threads)); // Fairly unbalanced algorithm
    apply->SetInput(0, ssfpData);
    apply->SetConst(0, T1);
    apply->SetConst(1, QI::ReadImage(args::get(B1)));
    apply->SetMask(QI::ReadImage(args::get(mask)));
    if (subregion) {
        apply->SetSubregion(QI::RegionOpt(args::get(subregion)));
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
    std::string outPrefix = args::get(outarg) + "FM_";
    QI::WriteImage(apply->GetOutput(0), outPrefix + "PD.nii");
    QI::WriteImage(apply->GetOutput(1), outPrefix + "T2.nii");
    QI::WriteImage(apply->GetOutput(2), outPrefix + "f0.nii");
    QI::WriteImage(apply->GetIterationsOutput(), outPrefix + "its.nii");
    QI::WriteScaledImage(apply->GetResidualOutput(), apply->GetOutput(0), outPrefix + "residual.nii");
    if (resids) {
        QI::WriteScaledVectorImage(apply->GetAllResidualsOutput(), apply->GetOutput(0), outPrefix + "all_residuals.nii");
    }
    return EXIT_SUCCESS;
}
