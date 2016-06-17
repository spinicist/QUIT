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

#include <getopt.h>
#include <iostream>

#include <Eigen/Dense>
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/NumericalDiff>

#include "cppoptlib/problem.h"
#include "cppoptlib/solver/lbfgsbsolver.h"
#include "Filters/ImageToVectorFilter.h"
#include "Filters/ApplyAlgorithmFilter.h"
#include "QI/Models/Models.h"
#include "QI/Sequences/Sequences.h"
#include "QI/Util.h"

using namespace std;
using namespace Eigen;

//******************************************************************************
// Algorithm Subclasses
//******************************************************************************
class D1Algo : public Algorithm<double> {
public:
    static const size_t DefaultIterations = 15;
protected:
    const shared_ptr<QI::Model> m_model = make_shared<QI::SCD>();
    shared_ptr<QI::SPGRSimple> m_sequence;
    size_t m_iterations = DefaultIterations;
    double m_thresh = -numeric_limits<double>::infinity();
    double m_lo = -numeric_limits<double>::infinity();
    double m_hi = numeric_limits<double>::infinity();

public:
    void setIterations(size_t n) { m_iterations = n; }
    size_t getIterations() { return m_iterations; }
    void setSequence(shared_ptr<QI::SPGRSimple> &s) { m_sequence = s; }
    void setThreshold(double t) { m_thresh = t; }
    void setClamp(double lo, double hi) { m_lo = lo; m_hi = hi; }
    size_t numInputs() const override { return m_sequence->count(); }
    size_t numConsts() const override { return 1; }
    size_t numOutputs() const override { return 2; }
    size_t dataSize() const override { return m_sequence->size(); }

    virtual TArray defaultConsts() override {
        // B1
        TArray def = TArray::Ones(1);
        return def;
    }
};

class D1LLS : public D1Algo {
public:
    virtual void apply(const TInput &data, const TArray &inputs,
                       TArray &outputs, TArray &resids, TIterations &its) const override
    {
        double B1 = inputs[0];
        ArrayXd flip = m_sequence->flip() * B1;
        VectorXd Y = data / flip.sin();
        MatrixXd X(Y.rows(), 2);
        X.col(0) = data / flip.tan();
        X.col(1).setOnes();
        VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
        outputs[1] = -m_sequence->TR() / log(b[0]);
        outputs[0] = b[1] / (1. - b[0]);
        ArrayXd theory = QI::One_SPGR(m_sequence->flip(), m_sequence->TR(), outputs[0], outputs[1], B1).array().abs();
        resids = (data.array() - theory);
        if (outputs[0] < m_thresh)
            outputs.setZero();
        outputs[1] = QI::clamp(outputs[1], m_lo, m_hi);
        its = 1;
    }
};

class D1WLLS : public D1Algo {
public:
    virtual void apply(const TInput &data, const TArray &inputs,
                       TArray &outputs, TArray &resids, TIterations &its) const override
    {
        double B1 = inputs[0];
        ArrayXd flip = m_sequence->flip() * B1;
        VectorXd Y = data / flip.sin();
        MatrixXd X(Y.rows(), 2);
        X.col(0) = data / flip.tan();
        X.col(1).setOnes();
        VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
        outputs[1] = -m_sequence->TR() / log(b[0]);
        outputs[0] = b[1] / (1. - b[0]);
        for (its = 0; its < m_iterations; its++) {
            VectorXd W = (flip.sin() / (1. - (exp(-m_sequence->TR()/outputs[1])*flip.cos()))).square();
            b = (X.transpose() * W.asDiagonal() * X).partialPivLu().solve(X.transpose() * W.asDiagonal() * Y);
            Array2d newOutputs;
            newOutputs[1] = -m_sequence->TR() / log(b[0]);
            newOutputs[0] = b[1] / (1. - b[0]);
            if (newOutputs.isApprox(outputs))
                break;
            else
                outputs = newOutputs;
        }
        ArrayXd theory = QI::One_SPGR(m_sequence->flip(), m_sequence->TR(), outputs[0], outputs[1], B1).array().abs();
        resids = (data.array() - theory);
        if (outputs[0] < m_thresh)
            outputs.setZero();
        outputs[1] = QI::clamp(outputs[1], m_lo, m_hi);
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
    using typename D1Algo::TArray;
    using typename D1Algo::TInput;
    using typename D1Algo::TIterations;

    virtual void apply(const TInput &indata, const TArray &consts, TArray &outputs, TArray &resids, TIterations &its) const override {
        double B1 = consts[0];
        // Improve scaling by dividing the PD down to something sensible.
        // This gets scaled back up at the end.
        const TInput data = indata / indata.maxCoeff();
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
        outputs[0] = p[0] * indata.maxCoeff();
        outputs[1] = p[1];
        if (!isfinite(p[1])) {
            cout << "Not finite" << endl;
            cout << indata.transpose() << endl;
            cout << data.transpose() << endl;
            cout << p.transpose() << endl << B1 << " " << indata.abs().maxCoeff() << endl;
            solver.setDebug(cppoptlib::DebugLevel::High);
            Eigen::Vector2d p; p << 10., 1.0;
            solver.minimize(cost, p);
        }
        resids = cost.residuals(p) * indata.abs().maxCoeff();
    }
};

class D1NLLS : public D1Algo {
public:
    virtual void apply(const TInput &indata, const TArray &inputs,
                       TArray &outputs, TArray &resids, TIterations &its) const override
    {
        double B1 = inputs[0];
        const TInput data = indata / indata.maxCoeff();
        T1Functor f(m_sequence, data, B1);
        NumericalDiff<T1Functor> nDiff(f);
        LevenbergMarquardt<NumericalDiff<T1Functor>> lm(nDiff);
        lm.setMaxfev(m_iterations * (m_sequence->size() + 1));
        // PD & T1 - Initial guess of 1s
        VectorXd p(2); p << data.maxCoeff() * 10., 1.;
        lm.minimize(p);
        outputs = p;
        ArrayXd theory = QI::One_SPGR(m_sequence->flip(), m_sequence->TR(), outputs[0], outputs[1], B1).array().abs();
        resids = indata.maxCoeff() * (data.array() - theory);
        outputs[0] *= indata.maxCoeff();
        if (outputs[0] < m_thresh)
            outputs.setZero();
        outputs[1] = QI::clamp(outputs[1], m_lo, m_hi);
        its = lm.iterations();
    }
};

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: qidespot1 [options] spgr_input \n\
\n\
Options:\n\
\n\
    --help, -h        : Print this message\n\
    --verbose, -v     : Print more information\n\
    --no-prompt, -n   : Suppress input prompts\n\
    --out, -o path    : Add a prefix to the output filenames\n\
    --mask, -m file   : Mask input with specified file\n\
    --B1, -b file     : B1 Map file (ratio)\n\
    --thresh, -t n    : Threshold maps at PD < n\n\
    --clamp, -c n     : Clamp T1 between 0 and n\n\
    --algo, -a l      : LLS algorithm (default)\n\
               w      : WLLS algorithm\n\
               n      : NLLS (Levenberg-Marquardt)\n\
               b      : LBFGSB algorithm\n\
    --its, -i N       : Max iterations for WLLS/NLLS (default 15)\n\
    --resids, -r      : Write out per flip-angle residuals\n\
    --threads, -T N   : Use N threads (default=4, 0=hardware limit)\n"
};

bool verbose = false, prompt = true, all_residuals = false;
int num_threads = 4;
string outPrefix;
const struct option long_options[] =
{
    {"help", no_argument, 0, 'h'},
    {"verbose", no_argument, 0, 'v'},
    {"no-prompt", no_argument, 0, 'n'},
    {"out", required_argument, 0, 'o'},
    {"mask", required_argument, 0, 'm'},
    {"B1", required_argument, 0, 'b'},
    {"thresh", required_argument, 0, 't'},
    {"clamp", required_argument, 0, 'c'},
    {"algo", required_argument, 0, 'a'},
    {"its", required_argument, 0, 'i'},
    {"resids", no_argument, 0, 'r'},
    {"threads", required_argument, 0, 'T'},
    {0, 0, 0, 0}
};
static const char *short_opts = "hvnm:o:b:t:c:a:i:rT:";

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();
    QI::VolumeF::Pointer mask = ITK_NULLPTR;
    QI::VolumeF::Pointer B1   = ITK_NULLPTR;

    shared_ptr<D1Algo> algo = make_shared<D1LLS>();
    int indexptr = 0, c;
    while ((c = getopt_long(argc, argv, short_opts, long_options, &indexptr)) != -1) {
        switch (c) {
            case 'v': verbose = true; break;
            case 'n': prompt = false; break;
            case 'a':
                switch (*optarg) {
                    case 'l': algo = make_shared<D1LLS>();  if (verbose) cout << "LLS algorithm selected." << endl; break;
                    case 'w': algo = make_shared<D1WLLS>(); if (verbose) cout << "WLLS algorithm selected." << endl; break;
                    case 'n': algo = make_shared<D1NLLS>(); if (verbose) cout << "NLLS algorithm selected." << endl; break;
                    case 'b': algo = make_shared<D1LBFGSB>(); if (verbose) cout << "LBFGSB algorithm selected." << endl; break;
                    default:
                        cerr << "Unknown algorithm type " << optarg << endl;
                        return EXIT_FAILURE;
                        break;
                } break;
            default: break;
        }
    }
    optind = 1;
    while ((c = getopt_long(argc, argv, short_opts, long_options, &indexptr)) != -1) {
        switch (c) {
            case 'v': case 'n': case 'a': break;
            case 'm':
                if (verbose) cout << "Opening mask file " << optarg << endl;
                mask = QI::ReadImage(optarg);
                break;
            case 'o':
                outPrefix = optarg;
                if (verbose) cout << "Output prefix will be: " << outPrefix << endl;
                break;
            case 'b':
                if (verbose) cout << "Opening B1 file: " << optarg << endl;
                B1 = QI::ReadImage(optarg);
                break;
            case 't': algo->setThreshold(atof(optarg)); break;
            case 'c': algo->setClamp(0, atof(optarg)); break;
            case 'i': algo->setIterations(atoi(optarg)); break;
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
    if ((argc - optind) != 1) {
        cout << "Incorrect number of arguments." << endl << usage << endl;
        return EXIT_FAILURE;
    }

    string inputFilename = argv[optind++];
    if (verbose) cout << "Opening SPGR file: " << inputFilename << endl;
    auto data = QI::ReadVectorImage<float>(inputFilename);
    shared_ptr<QI::SPGRSimple> spgrSequence = make_shared<QI::SPGRSimple>(cin, prompt);
    if (verbose) cout << *spgrSequence;
    algo->setSequence(spgrSequence);
    auto apply = itk::ApplyAlgorithmFilter<D1Algo>::New();
    apply->SetAlgorithm(algo);
    apply->SetOutputAllResiduals(all_residuals);
    apply->SetPoolsize(num_threads);
    apply->SetInput(0, data);
    if (mask)
        apply->SetMask(mask);
    if (B1)
        apply->SetConst(0, B1);
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
    outPrefix = outPrefix + "D1_";
    QI::WriteImage(apply->GetOutput(0), outPrefix + "PD.nii");
    QI::WriteImage(apply->GetOutput(1), outPrefix + "T1.nii");
    QI::WriteScaledImage(apply->GetResidualOutput(), apply->GetOutput(0), outPrefix + "residual.nii");
    if (all_residuals) {
        QI::WriteScaledVectorImage(apply->GetAllResidualsOutput(), apply->GetOutput(0), outPrefix + "all_residuals.nii");
    }
    if (algo->getIterations() != D1Algo::DefaultIterations) {
        QI::WriteImage(apply->GetIterationsOutput(), outPrefix + "iterations.nii");
    }
    if (verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}
