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
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/NumericalDiff>

#include "itkImageFileReader.h"

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

// HIFI Algorithm - includes optimising B1
class HIFIFunctor : public DenseFunctor<double> {
    protected:
        const shared_ptr<QI::SPGRSimple> m_spgr;
        const shared_ptr<QI::MPRAGE> m_mprage;
        const ArrayXd m_data;

    public:
        HIFIFunctor(const shared_ptr<QI::SPGRSimple> spgr,
                    const shared_ptr<QI::MPRAGE> mprage,
                    const ArrayXd &data) :
            DenseFunctor<double>(3, spgr->size() + mprage->size()),
            m_spgr(spgr), m_mprage(mprage), m_data(data)
        {
            assert(static_cast<size_t>(m_data.rows()) == values());
        }

        int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
            eigen_assert(diffs.size() == values());

            ArrayXd s(values());
            s.head(m_spgr->size()) = QI::One_SPGR_Magnitude(m_spgr->flip(), m_spgr->TR(), params[0], params[1], params[2]);
            s.tail(m_mprage->size()) = QI::One_MPRAGE(m_mprage->flip()[0], m_mprage->TR(), m_mprage->m_Nseg, m_mprage->m_Nk0, m_mprage->m_TI, m_mprage->m_TD,
                                                      params[0], params[1], params[2], m_mprage->m_eta).cwiseAbs(); 
            diffs = s - m_data;
            return 0;
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

    virtual void apply(const std::vector<TInput> &inputs, const std::vector<TConst> &, // No constants, remove name to silence compiler warnings
                       std::vector<TOutput> &outputs, TConst &residual,
                       TInput &resids, TIters &its) const override
    {
        Eigen::Map<const ArrayXf> spgr(inputs[0].GetDataPointer(), inputs[0].Size());
        Eigen::Map<const ArrayXf> irspgr(inputs[1].GetDataPointer(), inputs[1].Size());
        ArrayXf indata(dataSize()); indata << spgr, irspgr;
        const ArrayXd data = indata.cast<double>() / indata.abs().maxCoeff(); // Scale to make parameters roughly equal
        HIFIFunctor f(m_spgr, m_mprage, data);
        NumericalDiff<HIFIFunctor> nDiff(f);
        LevenbergMarquardt<NumericalDiff<HIFIFunctor>> lm(nDiff);
        // LevenbergMarquardt does not currently have a good interface, have to do things in steps
        lm.setMaxfev(m_iterations * (data.rows() + 1));
        VectorXd p(3); p << 10., 1., 1.; // Initial guess
        lm.minimize(p);
        ArrayXd r(indata.rows());
        f(p, r); // Get the residuals
        
        if (p[0] < m_thresh)
            p[0] = 0;
        else
            outputs[0] = p[0] * indata.abs().maxCoeff();
        outputs[1] = QI::clamp(p[1], m_lo, m_hi);
        outputs[2] = p[2];
        residual = sqrt(r.square().sum() / r.rows());
        ArrayXf rf = r.cast<float>();
        resids = itk::VariableLengthVector<float>(rf.data(), rf.rows());
        its = lm.iterations();
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
        cout << "Mean time per voxel was " << apply->GetMeanTime() << "s" << endl;
        cout << "Writing results files." << endl;
    }
    outPrefix = outPrefix + "HIFI_";

    QI::WriteImage(apply->GetOutput(0), outPrefix + "PD.nii");
    QI::WriteImage(apply->GetOutput(1), outPrefix + "T1.nii");
    QI::WriteImage(apply->GetOutput(2), outPrefix + "B1.nii");
    QI::WriteScaledImage(apply->GetResidualOutput(), apply->GetOutput(0), outPrefix + "residual.nii");
    if (all_residuals) {
        QI::WriteScaledVectorImage(apply->GetAllResidualsOutput(), apply->GetOutput(0), outPrefix + "all_residuals.nii");
    }
    if (verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}
