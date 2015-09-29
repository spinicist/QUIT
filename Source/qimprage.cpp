/*
 *  qimp2rage.cpp
 *
 *  Created by Tobias Wood on 2015/08/24.
 *  Copyright (c) 2015 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <string>
#include <complex>
#include <getopt.h>

#include "itkBinaryFunctorImageFilter.h"

#include "Filters/ApplyAlgorithmFilter.h"
#include "Types.h"
#include "Util.h"
#include "Sequence.h"

using namespace std;

const string usage{
"A tool to process MP3-RAGE images\n\
Usage is: qimp2rage [options] input \n\
\
Input must contain three volumes and be complex-valued\n\
Output will be files with _CXY where X & Y are inversion times\n\
Options:\n\
    --verbose, -v     : Print more messages.\n\
    --mask, -m file   : Mask input with specified file.\n\
    --out, -o path    : Add a prefix to the output filenames.\n\
    --complex, -x     : Output complex contast images.\n\
    --quant, -q       : Generate T1, B1 and inversion efficiency maps\n\
    --threads, -T N   : Use a maximum of N threads.\n"
};

static struct option long_options[] = {
    {"verbose", required_argument, 0, 'v'},
    {"mask",    required_argument, 0, 'm'},
    {"out",     required_argument, 0, 'o'},
    {"complex", no_argument, 0, 'x'},
    {"quant", no_argument, 0, 'q'},
    {"threads", required_argument, 0, 'T'},
    {0, 0, 0, 0}
};
static const char *short_options = "vm:o:xqT:h";

template<class T> class MP2Contrast {
public:
    MP2Contrast() {}
    ~MP2Contrast() {}
    bool operator!=(const MP2Contrast &) const { return false; }
    bool operator==(const MP2Contrast &other) const { return !(*this != other); }

    inline complex<T> operator()(const complex<T> &ti1, const complex<T> &ti2) const
    {
        const T a1 = abs(ti1);
        const T a2 = abs(ti2);
        return (conj(ti1)*ti2)/(a1*a1 + a2*a2);
    }
};

class MPRAGEAlgo : public Algorithm<complex<double>> {
protected:
    std::vector<Array3d> m_pars, m_cons;
    bool m_do_lookup = false;

public:
    size_t numInputs() const override { return 1; }
    size_t numConsts() const override { return 0; }
    size_t numOutputs() const override { return 6; }
    size_t dataSize() const override { return 3; }
    void setDoLookup(const bool l) {
        m_do_lookup = l;
        if (m_do_lookup) {
            MP3RAGE sequence(true);
            MP2Contrast<double> con;
            m_pars.clear();
            m_cons.clear();
            for (float T1 = 0.5; T1 < 2.5; T1 += 0.01) {
                for (float B1 = 0.75; B1 < 1.25; B1 += 0.05) {
                    for (float eta = 0.5; eta < 1.5; eta += 0.05) {
                        Array3d tp; tp << T1, B1, eta;
                        m_pars.push_back(tp);
                        Array3cd sig = sequence.signal(1., T1, B1, eta);
                        Array3d tc;
                        tc[0] = con(sig[0], sig[1]).real();
                        tc[1] = con(sig[0], sig[2]).real();
                        tc[2] = con(sig[1], sig[2]).real();
                        m_cons.push_back(tc);
                        //cout << m_pars.back().transpose() << " : " << m_cons.back().transpose() << endl;
                    }
                }
            }
            cout << "Lookup table has " << m_pars.size() << " entries" << endl;
        }
    }

    virtual TArray defaultConsts() override {
        // B1
        TArray def = TArray::Ones(0);
        return def;
    }

    MPRAGEAlgo() {
    }

    virtual void apply(const TInput &data_inputs, const TArray &const_inputs,
                       TArray &outputs, TArray &resids, TIterations &its) const override {
        double best_distance = numeric_limits<double>::max();
        int best_index = 0;

        MP2Contrast<double> con;
        Array3d in_cons;
        in_cons[0] = con(data_inputs[0], data_inputs[1]).real();
        in_cons[1] = con(data_inputs[0], data_inputs[2]).real();
        in_cons[2] = con(data_inputs[1], data_inputs[2]).real();
        outputs.head(3) = in_cons;

        if (m_do_lookup) {
            for (int i = 0; i < m_pars.size(); i++) {
                double distance = (m_cons[i] - in_cons).matrix().norm();
                if (distance < best_distance) {
                    best_distance = distance;
                    best_index = i;
                }
            }
            outputs.tail(3) = m_pars[best_index];
            resids = in_cons - m_cons[best_index];
        } else {
            outputs.tail(3).setConstant(0);
            resids.setConstant(0);
        }
        //cout << "Best index " << best_index << " distance " << best_distance << " pars " << outputs.transpose() << " data " << data_inputs.transpose() << " cons" << m_cons[best_index].transpose() << endl;
        its = 1;
    }
};

int main(int argc, char **argv) {
    int indexptr = 0, c;
    string outPrefix = "";
    bool verbose = false, complex_output = false, do_lookup = false;
    QI::ImageF::Pointer mask = ITK_NULLPTR;
    while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
        switch (c) {
            case 'v': verbose = true; break;
            case 'm': {
                cout << "Reading mask." << endl;
                auto maskFile = QI::ReadImageF::New();
                maskFile->SetFileName(optarg);
                maskFile->Update();
                mask = maskFile->GetOutput();
                mask->DisconnectPipeline();
            } break;
            case 'o':
                outPrefix = optarg;
                cout << "Output prefix will be: " << outPrefix << endl;
                break;
            case 'x': complex_output = true; break;
            case 'q': do_lookup = true; break;
            case 'T': itk::MultiThreader::SetGlobalDefaultNumberOfThreads(atoi(optarg)); break;
            case 'h':
                cout << usage << endl;
                return EXIT_SUCCESS;
            case '?': // getopt will print an error message
                return EXIT_FAILURE;
            default:
                cout << "Unhandled option " << string(1, c) << endl;
                return EXIT_FAILURE;
        }
    }
    if ((argc - optind) != 1) {
        cout << usage << endl;
        return EXIT_FAILURE;
    }
    string fname(argv[optind]);
    if (verbose) cout << "Opening input file " << fname << endl;
    auto inFile = QI::ReadTimeseriesXF::New();
    inFile->SetFileName(fname);
    inFile->Update();
    if (verbose) cout << "Processing" << endl;
    auto apply = itk::ApplyAlgorithmFilter<MPRAGEAlgo, complex<float>>::New();
    auto process = make_shared<MPRAGEAlgo>();
    process->setDoLookup(do_lookup);
    apply->SetAlgorithm(process);

    auto vectorFilter = QI::TimeseriesToVectorXF::New();
    vectorFilter->SetInput(inFile->GetOutput());
    apply->SetInput(0, vectorFilter->GetOutput());
    apply->SetMask(mask);
    apply->Update();
    QI::writeResult(apply->GetOutput(0), outPrefix + "MP3_C12" + QI::OutExt());
    QI::writeResult(apply->GetOutput(1), outPrefix + "MP3_C13" + QI::OutExt());
    QI::writeResult(apply->GetOutput(2), outPrefix + "MP3_C23" + QI::OutExt());
    if (do_lookup) {
        QI::writeResult(apply->GetOutput(3), outPrefix + "MP3_T1" + QI::OutExt());
        QI::writeResult(apply->GetOutput(4), outPrefix + "MP3_B1" + QI::OutExt());
        QI::writeResult(apply->GetOutput(5), outPrefix + "MP3_eta" + QI::OutExt());
        QI::writeResiduals(apply->GetResidOutput(), outPrefix + "MP3_", true);
    }
    if (verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}

