/*
 *  qmcdespot.cpp
 *
 *  Created by Tobias Wood on 03/06/2015.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <getopt.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/NumericalDiff>

#include "itkTimeProbe.h"

#include "QI/Types.h"
#include "QI/Util.h"
#include "QI/IO.h"
#include "QI/Models/Model.h"
#include "QI/Sequences/Sequences.h"
#include "QI/RegionContraction.h"
#include "Filters/ApplyAlgorithmFilter.h"
#include "Filters/ReorderVectorFilter.h"

using namespace std;
using namespace Eigen;

/*
 * Read in all required files and data from cin
 */
void parseInput(shared_ptr<QI::SequenceGroup> seq, vector<typename QI::VectorVolumeF::Pointer> &images,
                bool flip, bool verbose, bool prompt);
void parseInput(shared_ptr<QI::SequenceGroup> seq, vector<typename QI::VectorVolumeF::Pointer> &images,
                bool flip, bool verbose, bool prompt)
{
    string path;
    if (prompt) cout << "Enter input filename: " << flush;
    while (QI::Read(cin, path) && (path != "END") && (path != "")) {
        if (verbose) cout << "Reading file: " << path << endl;
        auto image = QI::ReadVectorImage<float>(path);
        seq->addSequence(QI::ReadSequence(cin, prompt));
        image->DisconnectPipeline(); // This step is really important.
        images.push_back(image);
        if (prompt) cout << "Enter next filename (END to finish input): " << flush;
    }
}

class MCDAlgo : public QI::ApplyF::Algorithm {
protected:
    ArrayXXd m_bounds;
    shared_ptr<QI::Model> m_model = nullptr;
    shared_ptr<QI::SequenceGroup> m_sequence = nullptr;
    QI::FieldStrength m_tesla = QI::FieldStrength::Three;
    int m_iterations = 0;

public:
    MCDAlgo(shared_ptr<QI::Model>&m, ArrayXXd &b,
            shared_ptr<QI::SequenceGroup> s, int mi) :
        m_model(m), m_bounds(b), m_sequence(s), m_iterations(mi)
    {}

    size_t numInputs() const override  { return m_sequence->count(); }
    size_t numOutputs() const override { return m_model->nParameters(); }
    size_t dataSize() const override   { return m_sequence->size(); }

    void setModel(shared_ptr<QI::Model> &m) { m_model = m; }
    void setSequence(shared_ptr<QI::SequenceGroup> &s) { m_sequence = s; }
    void setBounds(ArrayXXd &b) { m_bounds = b; }
    void setIterations(const int i) { m_iterations = i; }
    const float &zero(const size_t i) const override { static float zero = 0; return zero; }
};

class MCDSRCFunctor {
    public:
        const shared_ptr<QI::SequenceGroup> m_sequence;
        const ArrayXd m_data, m_weights;
        const shared_ptr<QI::Model> m_model;

        MCDSRCFunctor(shared_ptr<QI::Model> m,shared_ptr<QI::SequenceGroup> s, const ArrayXd &d, const ArrayXd &w) :
            m_sequence(s), m_data(d), m_model(m), m_weights(w)
        {
            assert(static_cast<size_t>(m_data.rows()) == m_sequence->size());
        }

        int inputs() const { return m_model->nParameters(); }
        int values() const { return m_sequence->size(); }

        const bool constraint(const VectorXd &params) const {
            return m_model->ValidParameters(params);
        }

        ArrayXd residuals(const Ref<VectorXd> &params) const {
            const ArrayXd s = (m_sequence->signal(m_model, params)).abs();
            return m_data - s;
        }
        double operator()(const Ref<VectorXd> &params) const {
            return (residuals(params) * m_weights).square().sum();
        }
};

class SRCAlgo : public MCDAlgo {
using MCDAlgo::MCDAlgo;

private:
    size_t m_samples = 5000, m_retain = 50;
    bool m_gauss = true;

public:
    void setGauss(bool g) { m_gauss = g; }

    size_t numConsts() const override  { return 2; }
    virtual std::vector<float> defaultConsts() const override {
        std::vector<float> def(2);
        def[0] = NAN; def[1] = 1.0f; // f0, B1
        return def;
    }

virtual bool apply(const std::vector<TInput> &inputs, const std::vector<TConst> &consts,
                   std::vector<TOutput> &outputs, TConst &residual,
                   TInput &resids, TIters &its) const override
    {
        ArrayXd data(dataSize());
        int dataIndex = 0;
        for (int i = 0; i < inputs.size(); i++) {
            Eigen::Map<const ArrayXf> this_data(inputs[i].GetDataPointer(), inputs[i].Size());
            if (m_model->scaleToMean()) {
                data.segment(dataIndex, this_data.rows()) = this_data.cast<double>() / this_data.abs().mean();
            } else {
                data.segment(dataIndex, this_data.rows()) = this_data.cast<double>();
            }
            dataIndex += this_data.rows();
        }
        ArrayXd thresh(m_model->nParameters()); thresh.setConstant(0.05);
        const double f0 = consts[0];
        const double B1 = consts[1];
        ArrayXXd localBounds = m_bounds;
        ArrayXd weights = ArrayXd::Ones(m_sequence->size());
        if (isfinite(f0)) { // We have an f0 map, add it to the fitting bounds
            localBounds.row(m_model->ParameterIndex("f0")) += f0;
            weights = m_sequence->weights(f0);
        }
        localBounds.row(m_model->ParameterIndex("B1")).setConstant(B1);
        MCDSRCFunctor func(m_model, m_sequence, data, weights);
        QI::RegionContraction<MCDSRCFunctor> rc(func, localBounds, thresh, m_samples, m_retain, m_iterations, 0.02, m_gauss, false);
        ArrayXd pars(m_model->nParameters());
        rc.optimise(pars);
        for (int i = 0; i < m_model->nParameters(); i++) {
            outputs[i] = pars[i];
        }
        ArrayXf r = func.residuals(pars).cast<float>();
        residual = sqrt(r.square().sum() / r.rows());
        resids = itk::VariableLengthVector<float>(r.data(), r.rows());
        its = rc.contractions();
        return true;
    }
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();

    auto tesla = QI::FieldStrength::Three;
    int start_slice = 0, stop_slice = 0;
    int max_its = 4, num_threads = 4;
    int verbose = false, prompt = true, all_residuals = false, flipData = false;
    string outPrefix;
    enum class Algos { SRC, GRC };
    Algos which_algo = Algos::GRC;

    QI::VolumeF::Pointer mask, B1, f0 = ITK_NULLPTR;
    shared_ptr<QI::Model> model = make_shared<QI::MCD3>();
    typedef itk::VectorImage<float, 2> VectorSliceF;
    auto apply = QI::ApplyF::New();
    const string usage {
"Usage is: qimcdespot [options]\n\
\n\
The program will prompt for input (unless --no-prompt specified)\n\
\n\
All times (TR) are in SECONDS. All angles are in degrees.\n\
\n\
Options:\n\
    --help, -h        : Print this message\n\
    --verbose, -v     : Print more information\n\
    --no-prompt, -n   : Don't print prompts for input\n\
    --mask, -m file   : Mask input with specified file\n\
    --out, -o path    : Add a prefix to the output filenames\n\
    --model, -M 1     : Use 1 component model\n\
                2     : Use 2 component model\n\
                2nex  : Use 2 component, no exchange model\n\
                3     : Use 3 component model (default)\n\
                3_f0  : Use 3 components with a myelin resonance frequency\n\
                3nex  : Use 3 component, no exchange model\n\
    --f0, -f file     : Use f0 Map file (in Hertz)\n\
    --B1, -b file     : B1 Map file (ratio)\n\
    --scale, -S       : Normalise signals to mean\n\
    --algo, -a S      : Use Uniform distribution for Region Contraction\n\
               G      : Use Gaussian distribution for RC (default)\n\
    --iters, -i N     : Specify maximum number of iterations (default 4)\n\
    --flip, -F        : Data order is phase then flip-angle (default opposite)\n\
    --tesla, -t 3     : Boundaries suitable for 3T (default)\n\
                7     : Boundaries suitable for 7T \n\
                u     : User specified boundaries from stdin\n\
    --resids, -r      : Write out per flip-angle residuals\n\
    --threads, -T N   : Use N threads (default=4, 0=hardware limit)\n\
    -s \"I J K SI SJ SK\" : Only process a subregion starting at voxel I,J,K\n\
                            with size SI,SJ,SK. Must fit within image.\n"
    };

    const struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"verbose", no_argument, 0, 'v'},
        {"mask", required_argument, 0, 'm'},
        {"out", required_argument, 0, 'o'},
        {"f0", required_argument, 0, 'f'},
        {"B1", required_argument, 0, 'b'},
        {"scale", no_argument, 0, 'S'},
        {"algo", required_argument, 0, 'a'},
        {"iterations", required_argument, 0, 'i'},
        {"flip", no_argument, 0, 'F'},
        {"tesla", required_argument, 0, 't'},
        {"resids", no_argument, 0, 'r'},
        {"threads", required_argument, 0, 'T'},
        {"no-prompt", no_argument, 0, 'n'},
        {"model", required_argument, 0, 'M'},
        {0, 0, 0, 0}
    };
    const char* short_options = "hvm:o:f:b:s:Sa:t:FT:rnM:i:j:";

    // Deal with these options in first pass to ensure the correct model is selected
    int indexptr = 0, c;
    while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
        switch (c) {
            case 'v': verbose = true; break;
            case 'n': prompt = false; break;
            case 'M': {
                string choose_model(optarg);
                if (choose_model == "1") { model = make_shared<QI::SCD>(); }
                else if (choose_model == "2") { model = make_shared<QI::MCD2>(); }
                else if (choose_model == "2nex") { model = make_shared<QI::MCD2_NoEx>(); }
                else if (choose_model == "3") { model = make_shared<QI::MCD3>(); }
                else if (choose_model == "3_f0") { model = make_shared<QI::MCD3_f0>(); }
                else if (choose_model == "3nex") { model = make_shared<QI::MCD3_NoEx>(); }
            } break;
            default:
                break;
        }
    }
    // Now reset and do a second pass
    optind = 1;
    while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
        switch (c) {
            case 'v': case 'n': case 'M': break; // Already handled
            case 'm':
                if (verbose) cout << "Reading mask file " << optarg << endl;
                mask = QI::ReadImage(optarg);
                break;
            case 'o':
                outPrefix = optarg;
                if (verbose) cout << "Output prefix will be: " << outPrefix << endl;
                break;
            case 'f':
                if (verbose) cout << "Reading f0 file: " << optarg << endl;
                f0 = QI::ReadImage(optarg);
                break;
            case 'b':
                if (verbose) cout << "Reading B1 file: " << optarg << endl;
                B1 = QI::ReadImage(optarg);
                break;
            case 's': {
                ArrayXd vals; QI::ReadArray(optarg, vals);
                if (vals.rows() != 6) {
                    QI_EXCEPTION( "Subregion must have 3 start indices and 3 sizes." );
                }
                QI::ApplyF::TRegion::IndexType start;
                QI::ApplyF::TRegion::SizeType size;
                start[0] = vals[0]; start[1] = vals[1]; start[2] = vals[2];
                size[0]  = vals[3]; size[1] =  vals[4]; size[2]  = vals[5];
                QI::ApplyF::TRegion subregion;
                subregion.SetIndex(start);
                subregion.SetSize(size);
                apply->SetSubregion(subregion);
            } break;
            case 'S':
                if (verbose) cout << "Mean scaling selected." << endl;
                model->setScaleToMean(true);
                break;
            case 'a':
                switch (*optarg) {
                case 'S': which_algo = Algos::SRC; break;
                case 'G': which_algo = Algos::GRC; break;
                default:
                    cerr << "Unknown algorithm type " << *optarg << endl;
                    return EXIT_FAILURE;
                    break;
                } break;
            case 'F': flipData = true; if (verbose) cout << "Data order is phase, then flip-angle" << endl; break;
            case 'T': 
                num_threads = stoi(optarg);
                if (num_threads == 0)
                    num_threads = std::thread::hardware_concurrency();
                break;
            case 't':
                switch (*optarg) {
                    case '3': tesla = QI::FieldStrength::Three; break;
                    case '7': tesla = QI::FieldStrength::Seven; break;
                    case 'u': tesla = QI::FieldStrength::User; break;
                    default:
                        cerr << "Unknown boundaries type " << *optarg << endl;
                        return EXIT_FAILURE;
                        break;
                } break;
            case 'i': max_its = atoi(optarg); break;
            case 'r': all_residuals = true; break;
            case 'h':
                cout << QI::GetVersion() << endl << usage << endl;
                return EXIT_SUCCESS;
            case '?': // getopt will print an error message
                return EXIT_FAILURE;
            case 0: break; // Just a flag
            default:
                cout << "Unhandled option " << string(1, c) << endl;
                return EXIT_FAILURE;
        }
    }
    if ((argc - optind) != 0) {
        cerr << usage << endl << "Incorrect number of arguments." << endl;
        return EXIT_FAILURE;
    } else if (prompt) {
        cout << "Starting qimcdespot" << endl;
        cout << "Run with -h switch to see usage" << endl;
    }
    Array2d f0Bandwidth;

    shared_ptr<QI::SequenceGroup> sequences = make_shared<QI::SequenceGroup>();
    // Build a Functor here so we can query number of parameters etc.
    if (verbose) cout << "Using " << model->Name() << " model." << endl;
    vector<QI::VectorVolumeF::Pointer> images;
    parseInput(sequences, images, flipData, verbose, prompt);

    ArrayXXd bounds = model->Bounds(tesla);
    ArrayXd start = model->Default(tesla);
    if (tesla == QI::FieldStrength::User) {
        ArrayXd temp;
        if (prompt) cout << "Enter lower bounds" << endl;
        QI::ReadArray(cin, temp);
        bounds.col(0) = temp;
        if (prompt) cout << "Enter upper bounds" << endl;
        QI::ReadArray(cin, temp);
        bounds.col(1) = temp;
    }
    switch (which_algo) {
    case Algos::SRC: {
        if (verbose) cout << "Using SRC algorithm" << endl;
        shared_ptr<SRCAlgo> algo = make_shared<SRCAlgo>(model, bounds, sequences, max_its);
        algo->setGauss(false);
        apply->SetAlgorithm(algo);
    } break;
    case Algos::GRC: {
        if (verbose) cout << "Using GRC algorithm" << endl;
        shared_ptr<SRCAlgo> algo = make_shared<SRCAlgo>(model, bounds, sequences, max_its);
        algo->setGauss(true);
        apply->SetAlgorithm(algo);
    } break;
    }
    apply->SetOutputAllResiduals(all_residuals);
    apply->SetVerbose(verbose);
    apply->SetPoolsize(num_threads);
    apply->SetSplitsPerThread(num_threads); // mcdespot with a mask & threads is a very unbalanced algorithm
    for (int i = 0; i < images.size(); i++) {
        apply->SetInput(i, images[i]);
    }
    if (f0) {
        f0->Update();
        apply->SetConst(0, f0);
    }
    if (B1) {
        B1->Update();
        apply->SetConst(1, B1);
    }
    if (mask) {
        mask->Update();
        apply->SetMask(mask);
    }

    // Need this here so the bounds.txt file will have the correct prefix
    outPrefix = outPrefix + model->Name() + "_";
    if (verbose) {
        cout << *sequences;
        cout << "Bounds:" << endl <<  bounds.transpose() << endl;
        ofstream boundsFile(outPrefix + "bounds.txt");
        boundsFile << "Names: ";
        for (size_t p = 0; p < model->nParameters(); p++) {
            boundsFile << model->ParameterNames()[p] << "\t";
        }
        boundsFile << endl << "Bounds: " << endl << bounds.transpose() << endl;
        boundsFile.close();
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
    for (int i = 0; i < model->nParameters(); i++) {
        QI::WriteImage(apply->GetOutput(i), outPrefix + model->ParameterNames()[i] + QI::OutExt());
    }
    QI::WriteScaledImage(apply->GetResidualOutput(), apply->GetOutput(0), outPrefix + "residual.nii");
    if (all_residuals) {
        QI::WriteScaledVectorImage(apply->GetAllResidualsOutput(), apply->GetOutput(0), outPrefix + "all_residuals.nii");
    }
    QI::WriteImage(apply->GetIterationsOutput(), outPrefix + "iterations" + QI::OutExt());
    return EXIT_SUCCESS;
}

