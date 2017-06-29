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

#include <getopt.h>
#include <iostream>
#include <Eigen/Dense>
#include "ceres/ceres.h"

#include "QI/Sequences/Sequences.h"
#include "QI/Util.h"
#include "QI/IO.h"
#include "Filters/ApplyAlgorithmFilter.h"

using namespace std;
using namespace Eigen;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: qidespot1hifi [options] spgr_input ir-spgr_input\n\
\
Options:\n\
    --help, -h        : Print this message\n\
    --verbose, -v     : Print more information\n\
    --no-prompt, -n   : Suppress input prompts\n\
    --mprage, -M      : Use a generic MP-RAGE sequence, not GE IR-SPGR\n\
    --out, -o path    : Add a prefix to the output filenames\n\
    --mask, -m file   : Mask input with specified file\n\
    --thresh, -t n    : Threshold maps at PD < n\n\
    --clamp, -c n     : Clamp T1 between 0 and n\n\
    --its, -i N       : Max iterations for NLLS (default 4)\n\
    --resids, -r      : Write out per flip-angle residuals\n\
    --threads, -T N   : Use N threads (default=4, 0=hardware limit)\n"
};

static bool verbose = false, prompt = true, IR = true, all_residuals = false;
static size_t nIterations = 4, num_threads = 4;
static string outPrefix;
static double thresh = -numeric_limits<double>::infinity();
static double clamp_lo = -numeric_limits<double>::infinity(), clamp_hi = numeric_limits<double>::infinity();
static const struct option long_opts[] = {
    {"help", no_argument, 0, 'h'},
    {"verbose", no_argument, 0, 'v'},
    {"no-prompt", no_argument, 0, 'n'},
    {"mprage", no_argument, 0, 'M'},
    {"mask", required_argument, 0, 'm'},
    {"out", required_argument, 0, 'o'},
    {"thresh", required_argument, 0, 't'},
    {"clamp", required_argument, 0, 'c'},
    {"its", required_argument, 0, 'i'},
    {"resids", no_argument, 0, 'r'},
    {"threads", required_argument, 0, 'T'},
    {0, 0, 0, 0}
};
static const char *short_opts = "hvnMm:o:t:c:s:p:i:rT:";

class SPGRCost : public ceres::CostFunction {
protected:
    const QI::SPGRSimple &m_seq;
    const ArrayXd m_data;

public:
    SPGRCost(const QI::SPGRSimple &s, const ArrayXd &data) :
        m_seq(s), m_data(data)
    {
        mutable_parameter_block_sizes()->push_back(3);
        set_num_residuals(data.size());
    }

    virtual bool Evaluate(double const* const* p,
                          double* resids,
                          double** jacobians) const override
    {
        const double &M0 = p[0][0];
        const double &T1 = p[0][1];
        const double &B1 = p[0][2];

        const ArrayXd sa = sin(B1 * m_seq.m_flip);
        const ArrayXd ca = cos(B1 * m_seq.m_flip);
        const double E1 = exp(-m_seq.m_TR / T1);
        const ArrayXd denom = (1.-E1*ca);
        
        Eigen::Map<Eigen::ArrayXd> r(resids, m_data.size());
        r = M0*sa*(1-E1)/denom - m_data;
        
        // std::cout << "SPGR RESIDS" << std::endl;
        // std::cout << r.transpose() << std::endl;
        if (jacobians && jacobians[0]) {
            Eigen::Map<Eigen::Matrix<double, -1, -1, RowMajor>> j(jacobians[0], m_data.size(), 3);
            j.col(0) = (1-E1)*sa/denom;
            j.col(1) = E1*M0*m_seq.m_TR*(ca-1.)*sa/((denom*T1).square());
            j.col(2) = M0*m_seq.m_flip*(1.-E1)*(ca-E1)/denom.square();
        }
        return true;
    }
};

// Use AutoDiff for this
class IRCostFunction  {
protected:
    const QI::MPRAGE &m_seq;
    const ArrayXd m_data;

public:
    IRCostFunction(const QI::MPRAGE &s, const ArrayXd &data) :
        m_seq(s), m_data(data)
    {
    }

    template<typename T> bool operator() (const T* const p1, T* r) const
    {
        const T &M0 = p1[0];
        const T &T1 = p1[1];
        const T &B1 = p1[2];
        const double eta = -1.0; // Inversion efficiency defined as -1 < eta < 0

        const double TIs = m_seq.m_TI[0] - m_seq.m_TR*m_seq.m_Nk0; // Adjust TI for k0
        const T T1s = 1. / (1./T1 - log(cos(m_seq.m_flip[0] * B1))/m_seq.m_TR);
        const T M0s = M0 * (1. - exp(-m_seq.m_TR/T1)) / (1. - exp(-m_seq.m_TR/T1s));
        const T A_1 = M0s*(1. - exp(-(m_seq.m_Nseg*m_seq.m_TR)/T1s));

        const T A_2 = M0*(1. - exp(-m_seq.m_TD[0]/T1));
        const T A_3 = M0*(1. - exp(-TIs/T1));
        const T B_1 = exp(-(m_seq.m_Nseg*m_seq.m_TR)/T1s);
        const T B_2 = exp(-m_seq.m_TD[0]/T1);
        const T B_3 = eta*exp(-TIs/T1);

        const T A = A_3 + A_2*B_3 + A_1*B_2*B_3;
        const T B = B_1*B_2*B_3;
        const T M1 = A / (1. - B);

        r[0] = m_data[0] - (M0s + (M1 - M0s)*exp(-(m_seq.m_Nk0*m_seq.m_TR)/T1s)) * sin(m_seq.m_flip[0] * B1);
        // std::cout << "MPRAGE " << m_data[0] << " " << M0*(Ms + (M1 - Ms)*E1k) * sin(m_seq.m_flip[0] * b) << std::endl;
        return true;
    }
};

class HIFIAlgo : public QI::ApplyF::Algorithm {
private:
    shared_ptr<QI::SPGRSimple> m_spgr;
    shared_ptr<QI::MPRAGE> m_mprage;
    size_t m_iterations = 15; // From tests this seems to be a sensible maximum number
    double m_thresh = -numeric_limits<double>::infinity();
    double m_lo = -numeric_limits<double>::infinity();
    double m_hi = numeric_limits<double>::infinity();
public:
    void setSequences(const shared_ptr<QI::SPGRSimple> &s, const shared_ptr<QI::MPRAGE> &m) { m_spgr = s; m_mprage = m;}
    void setIterations(size_t n) { m_iterations = n; }
    void setThreshold(double t) { m_thresh = t; }
    void setClamp(double lo, double hi) { m_lo = lo; m_hi = hi; }
    size_t numInputs() const override  { return 2; }
    size_t numConsts() const override  { return 0; }
    size_t numOutputs() const override { return 3; }
    size_t dataSize() const override   { return m_spgr->size() + m_mprage->size(); }
    const float &zero(const size_t i) const override { static float zero = 0; return zero; }

    virtual std::vector<float> defaultConsts() const override {
        // No constants for HIFI
        std::vector<float> def(0);
        return def;
    }

    virtual bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &, // No constants, remove name to silence compiler warnings
                       std::vector<TOutput> &outputs, TConst &residual,
                       TInput &resids, TIters &its) const override
    {
        Eigen::Map<const Eigen::ArrayXf> spgr_in(inputs[0].GetDataPointer(), inputs[0].Size());
        Eigen::Map<const Eigen::ArrayXf> ir_in(inputs[1].GetDataPointer(), inputs[1].Size());
        double scale = std::max(spgr_in.maxCoeff(), ir_in.maxCoeff());
        const ArrayXd spgr_data = spgr_in.cast<double>() / scale;
        const ArrayXd ir_data = ir_in.cast<double>() / scale;
        double spgr_pars[] = {10., 1., 1.}; // PD, T1, B1
        ceres::Problem problem;
        problem.AddResidualBlock(new SPGRCost(*m_spgr, spgr_data), NULL, spgr_pars);
        ceres::CostFunction *IRCost = new ceres::AutoDiffCostFunction<IRCostFunction, 1, 3>(new IRCostFunction(*m_mprage, ir_data));
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
        outputs[1] = spgr_pars[1];
        outputs[2] = spgr_pars[2];
        if (!summary.IsSolutionUsable()) {
            std::cout << summary.FullReport() << std::endl;
        }
        its = summary.iterations.size();
        residual = summary.final_cost * scale;
        if (resids.Size() > 0) {
            assert(resids.Size() == data.size());
            vector<double> r_temp(spgr_data.size() + 1);
            problem.Evaluate(ceres::Problem::EvaluateOptions(), NULL, &r_temp, NULL, NULL);
            for (int i = 0; i < r_temp.size(); i++) {
                resids[i] = r_temp[i];
            }
        }
        return true;
    }
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();
    QI::VolumeF::Pointer mask = ITK_NULLPTR;
    auto hifi = make_shared<HIFIAlgo>();
    int indexptr = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
        switch (c) {
            case 'v': verbose = true; break;
            case 'n': prompt = false; break;
            case 'M': IR = false; break;
            case 'm':
                if (verbose) cout << "Opening mask file: " << optarg << endl;
                mask = QI::ReadImage(optarg);
                break;
            case 'o':
                outPrefix = optarg;
                if (verbose) cout << "Output prefix will be: " << outPrefix << endl;
                break;
            case 't': hifi->setThreshold(atof(optarg)); break;
            case 'c': hifi->setClamp(0, atof(optarg)); break;
            case 'i': hifi->setIterations(atoi(optarg)); break;
            case 'r': all_residuals = true; break;
            case 'T':
                num_threads = stoi(optarg);
                if (num_threads == 0)
                    num_threads = std::thread::hardware_concurrency();
                break;
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
        cerr << "Incorrect number of arguments." << endl;
        cout << QI::GetVersion() << endl << usage << endl;
        return EXIT_FAILURE;
    }
    
    if (verbose) cout << "Opening SPGR file: " << argv[optind] << endl;
    auto spgrImg = QI::ReadVectorImage<float>(argv[optind++]);
    auto spgrSequence = make_shared<QI::SPGRSimple>(cin, prompt);
    if (verbose) cout << "Opening IR-SPGR file: " << argv[optind] << endl;
    auto irImg = QI::ReadVectorImage<float>(argv[optind++]);
    shared_ptr<QI::MPRAGE> irSequence;
    if (IR) {
        irSequence = make_shared<QI::IRSPGR>(cin, prompt);
    } else {
        irSequence = make_shared<QI::MPRAGE>(cin, prompt);
    }
    if (verbose) cout << *spgrSequence << endl << *irSequence << endl;
    auto apply = QI::ApplyF::New();
    hifi->setSequences(spgrSequence, irSequence);
    apply->SetAlgorithm(hifi);
    apply->SetOutputAllResiduals(all_residuals);
    apply->SetPoolsize(num_threads);
    apply->SetSplitsPerThread(num_threads);
    apply->SetVerbose(verbose);
    apply->SetInput(0, spgrImg);
    apply->SetInput(1, irImg);
    if (mask)
        apply->SetMask(mask);
    if (verbose) {
        cout << "Processing..." << endl;
        auto monitor = QI::GenericMonitor::New();
        apply->AddObserver(itk::ProgressEvent(), monitor);
    }
    apply->Update();
    if (verbose) {
        cout << "Elapsed time was " << apply->GetTotalTime() << "s" << endl;
        cout << "Writing results files." << endl;
    }
    outPrefix = outPrefix + "HIFI_";

    QI::WriteImage(apply->GetOutput(0), outPrefix + "PD" + QI::OutExt());
    QI::WriteImage(apply->GetOutput(1), outPrefix + "T1" + QI::OutExt());
    QI::WriteImage(apply->GetOutput(2), outPrefix + "B1" + QI::OutExt());
    QI::WriteScaledImage(apply->GetResidualOutput(), apply->GetOutput(0), outPrefix + "residual"  + QI::OutExt());
    if (all_residuals) {
        QI::WriteVectorImage(apply->GetAllResidualsOutput(), outPrefix + "all_residuals" + QI::OutExt());
    }
    if (verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}
