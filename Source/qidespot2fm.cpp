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

#include <getopt.h>
#include <iostream>
 
#include <Eigen/Dense>

#include "cppoptlib/meta.h"
#include "cppoptlib/problem.h"
#include "cppoptlib/solver/lbfgsbsolver.h"
#include "cppoptlib/solver/cmaessolver.h"
#include "QI/Util.h"
#include "QI/Models/Models.h"
#include "QI/Sequences/Sequences.h"
#include "Filters/ApplyAlgorithmFilter.h"

using namespace std;
using namespace Eigen;

class FMAlgo : public Algorithm<double> {
protected:
	shared_ptr<QI::SSFPSimple> m_sequence;
    bool m_symmetric = true;

public:
    typedef typename Algorithm<double>::TArray TArray;
    typedef typename Algorithm<double>::TInput TInput;
    typedef typename Algorithm<double>::TIterations TIterations;

    void setSequence(shared_ptr<QI::SSFPSimple> s) { m_sequence = s; }
    void setSymmetric(const bool b) { m_symmetric = b; }
    
    size_t numInputs() const override  { return m_sequence->count(); }
    size_t numConsts() const override  { return 2; }
    size_t numOutputs() const override { return 3; }
    size_t dataSize() const override   { return m_sequence->size(); }

    virtual TArray defaultConsts() override {
        // T1 & B1
        TArray def = TArray::Ones(2);
        return def;
    }
};

class FMCostFunction : public cppoptlib::Problem<double> {
public:
    Eigen::ArrayXd m_data;
    double m_T1, m_B1;
    shared_ptr<QI::SSFPSimple> m_sequence;
    const shared_ptr<QI::SCD> m_model = make_shared<QI::SCD>();

    Eigen::ArrayXd residuals(const Eigen::VectorXd &p) {
        ArrayXd s = QI::One_SSFP_Echo_Magnitude(m_sequence->flip(), m_sequence->phase_incs(), m_sequence->TR(), p[0], m_T1, p[1], p[2], m_B1);
        Eigen::ArrayXd diff = s - m_data;
        return diff;
    }

    double value(const cppoptlib::Vector<double> &p) {
        return residuals(p).square().sum();;
    }
    
    void gradient(const  cppoptlib::Vector<double> &p,  cppoptlib::Vector<double> &grad) {
        ArrayXd  s   = QI::One_SSFP_Echo_Magnitude(m_sequence->flip(), m_sequence->phase_incs(), m_sequence->TR(), p[0], m_T1, p[1], p[2], m_B1);
        ArrayXXd drv = QI::One_SSFP_Echo_Derivs(m_sequence->flip(), m_sequence->phase_incs(), m_sequence->TR(), p[0], m_T1, p[1], p[2], m_B1);
        grad = 2*(drv.colwise()*(s - m_data)).colwise().sum();
    }
};

class BFGSAlgo : public FMAlgo {
public:
    using typename FMAlgo::TArray;
    using typename FMAlgo::TInput;
    using typename FMAlgo::TIterations;

    virtual void apply(const TInput &indata, const TArray &consts, TArray &outputs, TArray &resids, TIterations &its) const override {
        double T1 = consts[0];
        double B1 = consts[1];
        if (isfinite(T1) && (T1 > 0.001)) {
            // Improve scaling by dividing the PD down to something sensible.
            // This gets scaled back up at the end.
            const auto data = indata / indata.abs().maxCoeff();
            
            cppoptlib::LbfgsbSolver<double>::Info opts;
            opts.rate = 1.e-5;
            opts.iterations = 100;
            opts.gradNorm = 1.e-3;
            opts.m = 5;
            
            cppoptlib::LbfgsbSolver<double> solver(opts);
            FMCostFunction cost;
            cost.m_B1 = B1;
            cost.m_data = data;
            cost.m_sequence = this->m_sequence;
            cost.m_T1 = T1;
            Array3d lower; lower << 1.e-3, this->m_sequence->TR(), -0.75 / this->m_sequence->TR();
            Array3d upper; upper << 1.e3,  T1,                      0.75 / this->m_sequence->TR();
            if (this->m_symmetric)
                lower[2] = 0.;
            cost.setLowerBound(lower);
            cost.setUpperBound(upper);
            
            vector<double> f0_starts;
            if (this->m_symmetric) {
                f0_starts.push_back(5.);
                f0_starts.push_back(0.25 / this->m_sequence->TR());
            } else {
				f0_starts.push_back(0.);
                f0_starts.push_back( 0.25 / this->m_sequence->TR());
                f0_starts.push_back(-0.25 / this->m_sequence->TR());
                f0_starts.push_back( 0.5 / this->m_sequence->TR());
				f0_starts.push_back(-0.5 / this->m_sequence->TR());
            }
            double best = numeric_limits<double>::infinity();
            Eigen::Array3d bestP;
            its = 0;
            for (const double &f0 : f0_starts) {
                Eigen::VectorXd p(3); p << 10., 0.1 * T1, f0; // Yarnykh gives T2 = 0.045 * T1 in brain, but best to overestimate for CSF
                solver.minimize(cost, p);
                double r = cost(p);
                if (r < best) {
                    best = r;
                    bestP = p;
					its = solver.info().iterations;
                }
            }
            outputs[0] = bestP[0] * indata.abs().maxCoeff();
            outputs[1] = bestP[1];
            outputs[2] = bestP[2];
            resids = cost.residuals(bestP) * indata.abs().maxCoeff();
        } else {
            outputs.setZero();
            resids.setZero();
            its = 0;
        }
    }
};

class CMAlgo : public FMAlgo {
public:
    using typename FMAlgo::TArray;
    using typename FMAlgo::TInput;
    using typename FMAlgo::TIterations;

    virtual void apply(const TInput &indata, const TArray &consts, TArray &outputs, TArray &resids, TIterations &its) const override {
        double T1 = consts[0];
        double B1 = consts[1];
        if (isfinite(T1) && (T1 > 0.001)) {
            // Improve scaling by dividing the PD down to something sensible.
            // This gets scaled back up at the end.
            const auto data = indata / indata.abs().maxCoeff();
            
            cppoptlib::LbfgsbSolver<double>::Info opts;
            opts.rate = 1.e-5;
            opts.iterations = 50;
            opts.gradNorm = 1.e-3;
            opts.m = 5;
            
            cppoptlib::CMAesSolver<double> solver(opts);
            FMCostFunction cost;
            cost.m_B1 = B1;
            cost.m_data = data;
            cost.m_sequence = this->m_sequence;
            cost.m_T1 = T1;
            Array3d lower; lower << 0.001, this->m_sequence->TR(), -0.6 / this->m_sequence->TR();
            Array3d upper; upper << 1.e3,  T1,                      0.6 / this->m_sequence->TR();
            if (this->m_symmetric)
                lower[2] = 0.;
            cost.setLowerBound(lower);
            cost.setUpperBound(upper);
            
            vector<double> f0_starts;
            if (this->m_symmetric) {
                f0_starts.push_back(5.);
                f0_starts.push_back(0.25 / this->m_sequence->TR());
            } else {
                f0_starts.push_back(0.);
                f0_starts.push_back(0.25 / this->m_sequence->TR());
                f0_starts.push_back(-0.25/ this->m_sequence->TR());
            }
            double best = numeric_limits<double>::infinity();
            Eigen::Array3d bestP;
            its = 0;
            for (const double &f0 : f0_starts) {
                Eigen::VectorXd p(3); p << 10., 0.1 * T1, f0; // Yarnykh gives T2 = 0.045 * T1 in brain, but best to overestimate for CSF
                solver.minimize(cost, p);
                double r = cost(p);
                if (r < best) {
                    best = r;
                    bestP = p;
                }
                its += solver.info().iterations;
            }
            outputs[0] = bestP[0] * indata.abs().maxCoeff();
            outputs[1] = bestP[1];
            outputs[2] = bestP[2];
            resids = cost.residuals(bestP) * indata.abs().maxCoeff();
        } else {
            outputs.setZero();
            resids.setZero();
            its = 0;
        }
    }
};


const string usage {
"Usage is: qidespot2fm [options] T1_map ssfp_file\n\
\
Options:\n\
    --help, -h        : Print this message\n\
    --verbose, -v     : Print more information\n\
    --no-prompt, -n   : Suppress input prompts\n\
    --mask, -m file   : Mask input with specified file\n\
    --out, -o path    : Add a prefix to the output filenames\n\
    --B1, -b file     : B1 Map file (ratio)\n\
    --algo, -a b      : Use LBFGSB algorithm (default)\n\
               c      : Use Covariance Maximization algorithm\n\
    --asym, -A        : Fit +/- off-resonance frequency\n\
    --flex, -f        : Specify all phase-incs for all flip-angles\n\
    --start, -s N     : Start processing from slice N\n\
    --stop, -p  N     : Stop processing at slice N\n\
    --finite, -F      : Use finite pulse length correction\n\
    --resids, -r      : Write out per flip-angle residuals\n\
    --threads, -T N   : Use N threads (default=4, 0=hardware limit)\n"
};

struct option long_opts[] = {
    {"help", no_argument, 0, 'h'},
    {"verbose", no_argument, 0, 'v'},
    {"no-prompt", no_argument, 0, 'n'},
    {"mask", required_argument, 0, 'm'},
    {"out", required_argument, 0, 'o'},
    {"B1", required_argument, 0, 'b'},
    {"algo", required_argument, 0, 'a'},
    {"asym", no_argument, 0, 'A'},
    {"flex", no_argument, 0, 'f'},
    {"start", required_argument, 0, 's'},
    {"stop", required_argument, 0, 'p'},
    {"threads", required_argument, 0, 'T'},
    {"resids", no_argument, 0, 'r'},
    {0, 0, 0, 0}
};
const char* short_opts = "hvnm:o:b:a:Afs:p:T:rd:";
int indexptr = 0;
char c;

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	Eigen::initParallel();
    typedef itk::ApplyAlgorithmFilter<FMAlgo> TApply;

    int start_slice = 0, stop_slice = 0;
    int verbose = false, prompt = true, all_residuals = false, symmetric = true,
        fitFinite = false, flex = false, use_BFGS = true, num_threads = 4;
    bool useCM = false;
    string outPrefix;
    QI::VolumeF::Pointer mask = ITK_NULLPTR, B1 = ITK_NULLPTR;

    optind = 1;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
        switch (c) {
        case 'v': verbose = true; break;
        case 'n': prompt = false; break;
        case 'A': symmetric = false; if (verbose) cout << "Asymmetric f0 fitting selected." << endl; break;
        case 'a':
            switch (*optarg) {
                case 'b': useCM = false; if (verbose) cout << "LBFGSB algorithm selected." << endl; break;
                case 'c': useCM = true;  if (verbose) cout << "CovMax algorithm selected." << endl; break;
            } break;
        case 'm':
            if (verbose) cout << "Reading mask file " << optarg << endl;
            mask = QI::ReadImage(optarg);
            break;
        case 'o':
            outPrefix = optarg;
            if (verbose) cout << "Output prefix will be: " << outPrefix << endl;
            break;
        case 'b':
            if (verbose) cout << "Reading B1 file: " << optarg << endl;
            B1 = QI::ReadImage(optarg);
            break;
        case 's': start_slice = atoi(optarg); break;
        case 'p': stop_slice = atoi(optarg); break;
        case 'f': flex = true; if (verbose) cout << "Flexible sequence input selected" << endl; break;
        case 'T':
            num_threads = stoi(optarg);
            if (num_threads == 0)
                num_threads = std::thread::hardware_concurrency();
            break;
        case 'r': all_residuals = true; break;
        case 0: break; // Just a flag
        case 'h':
            cout << QI::GetVersion() << endl << usage << endl;
            return EXIT_SUCCESS;
        case '?': // getopt will print an error message
            return EXIT_FAILURE;
        default:
            cout << "Unhandled option " << string(1, c) << endl;
            return EXIT_FAILURE;
        }
    }
    if ((argc - optind) != 2) {
        cout << QI::GetVersion() << endl << usage << endl;
        cout << "Wrong number of arguments. Need a T1 map and one SSFP file." << endl;
        return EXIT_FAILURE;
    }

    shared_ptr<QI::SSFPSimple> ssfpSequence;
    if (fitFinite) {
        cout << "Using finite pulse model." << endl;
        ssfpSequence = make_shared<QI::SSFPFinite>(cin, prompt);
    } else {
        if (flex)
            ssfpSequence = make_shared<QI::SSFPEchoFlex>(cin, prompt);
        else
            ssfpSequence = make_shared<QI::SSFPEcho>(cin, prompt);
    }
    if (verbose) cout << *ssfpSequence << endl;

    if (verbose) cout << "Reading T1 Map from: " << argv[optind] << endl;
    auto T1 = QI::ReadImage(argv[optind++]);
    if (verbose) cout << "Opening SSFP file: " << argv[optind] << endl;
    auto ssfpData = QI::ReadVectorImage<T>(argv[optind++]);
    auto apply = TApply::New();
    shared_ptr<FMAlgo> algo;
    if (useCM) {
        algo = make_shared<CMAlgo>();
    } else {
        algo = make_shared<BFGSAlgo>();
    }
    algo->setSequence(ssfpSequence);
    algo->setSymmetric(symmetric);
    apply->SetVerbose(verbose);
    apply->SetAlgorithm(algo);
    apply->SetPoolsize(num_threads);
    apply->SetInput(0, ssfpData);
    apply->SetConst(0, T1);
    apply->SetSlices(start_slice, stop_slice);
    if (B1) {
        apply->SetConst(1, B1);
    }
    if (mask) {
        apply->SetMask(mask);
    }
    if (verbose) {
        cout << "Processing" << endl;
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    if (verbose) {
        cout << "Elapsed time was " << apply->GetTotalTime() << "s" << endl;
        cout << "Mean time per voxel was " << apply->GetMeanTime() << "s" << endl;
        cout << "Writing results files." << endl;
    }
    outPrefix = outPrefix + "FM_";
    QI::WriteImage(apply->GetOutput(0), outPrefix + "PD.nii");
    QI::WriteImage(apply->GetOutput(1), outPrefix + "T2.nii");
    QI::WriteImage(apply->GetOutput(2), outPrefix + "f0.nii");
	QI::WriteImage(apply->GetIterationsOutput(), outPrefix + "its.nii");
    QI::WriteResiduals(apply->GetResidOutput(), outPrefix, all_residuals, apply->GetOutput(0));
    return EXIT_SUCCESS;
}
