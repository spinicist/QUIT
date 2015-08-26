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
#include "itkImageToHistogramFilter.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkComplexToRealImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkMaskImageFilter.h"

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
    --threshold, -t N : Threshold at Nth quantile.\n\
    --complex, -x     : Output complex contast images.\n\
    --threads, -T N   : Use a maximum of N threads.\n"
};

static struct option long_options[] = {
    {"verbose", required_argument, 0, 'v'},
    {"mask",    required_argument, 0, 'm'},
    {"threshold", required_argument, 0, 't'},
    {"out",     required_argument, 0, 'o'},
    {"complex", no_argument, 0, 'x'},
    {"threads", required_argument, 0, 'T'},
    {0, 0, 0, 0}
};
static const char *short_options = "vm:o:t:xT:h";

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

class MP3LookupAlgo : public Algorithm<double> {
protected:
    std::vector<Array3d> m_pars, m_cons;

public:
    size_t numInputs() const override { return 1; }
    size_t numConsts() const override { return 0; }
    size_t numOutputs() const override { return 3; }
    size_t dataSize() const override { return 3; }

    virtual TArray defaultConsts() {
        // B1
        TArray def = TArray::Ones(0);
        return def;
    }

    MP3LookupAlgo() {
        cout << "Build lookup table" << endl;
        MP3RAGE sequence(true);
        MP2Contrast<double> con;
        m_pars.clear();
        m_cons.clear();
        for (float T1 = 0.5; T1 < 4.3; T1 += 0.005) {
            for (float B1 = 0.75; B1 < 1.25; B1 += 0.01) {
                for (float eta = 1.0; eta < 1.25; eta += 1.0) {
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

    virtual void apply(const TInput &data_inputs, const TArray &const_inputs,
                       TArray &outputs, TArray &resids, TIterations &its) const override {
        double best_distance = numeric_limits<double>::max();
        int best_index = 0;

        for (int i = 0; i < m_pars.size(); i++) {
            double distance = (m_cons[i] - data_inputs).matrix().norm();
            if (distance < best_distance) {
                best_distance = distance;
                best_index = i;
            }
        }
        outputs = m_pars[best_index];
        //cout << "Best index " << best_index << " distance " << best_distance << " pars " << outputs.transpose() << " data " << data_inputs.transpose() << " cons" << m_cons[best_index].transpose() << endl;
        its = 1;
    }
};

int main(int argc, char **argv) {
    int indexptr = 0, c;
    string outPrefix = "";
    bool verbose = false, complex_output = false;
    float thresh_quantile = 0.0;
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
            case 't': thresh_quantile = atof(optarg); break;
            case 'x': complex_output = true; break;
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
    if (outPrefix == "")
        outPrefix = fname.substr(0, fname.find(".nii"));
    auto inFile = QI::ReadTimeseriesXF::New();
    inFile->SetFileName(fname);
    inFile->Update();
    vector<itk::ExtractImageFilter<QI::TimeseriesXF, QI::ImageXF>::Pointer> vols(3);

    auto region = inFile->GetOutput()->GetLargestPossibleRegion();
    region.GetModifiableSize()[3] = 0;

    for (int i = 0; i < 3; i++) {
        region.GetModifiableIndex()[3] = i;
        vols[i] = itk::ExtractImageFilter<QI::TimeseriesXF, QI::ImageXF>::New();
        vols[i]->SetExtractionRegion(region);
        vols[i]->SetInput(inFile->GetOutput());
        vols[i]->SetDirectionCollapseToSubmatrix();
    }

    if (!mask) {
        // Threshold the last volume to automatically generate a mask
        auto magFilter = itk::ComplexToModulusImageFilter<QI::ImageXF, QI::ImageF>::New();
        auto histFilter = itk::Statistics::ImageToHistogramFilter<QI::ImageF>::New();
        auto threshFilter = itk::BinaryThresholdImageFilter<QI::ImageF, QI::ImageF>::New();
        magFilter->SetInput(vols[2]->GetOutput());
        histFilter->SetInput(magFilter->GetOutput());
        itk::Statistics::ImageToHistogramFilter<QI::ImageF>::HistogramSizeType size(1); size.Fill(100);
        histFilter->SetHistogramSize(size);
        histFilter->SetAutoMinimumMaximum(true);
        histFilter->Update();
        float threshold = histFilter->GetOutput()->Quantile(0, thresh_quantile);
        threshFilter->SetInput(magFilter->GetOutput());
        threshFilter->SetLowerThreshold(threshold);
        threshFilter->SetInsideValue(1);
        threshFilter->SetOutsideValue(0);
        threshFilter->Update();
        mask = threshFilter->GetOutput();
        mask->DisconnectPipeline();
    }

    cout << "Generating contrast images" << endl;
    vector<QI::ImageF::Pointer> conImages(3, ITK_NULLPTR);
    int ind = 0;
    for (int i1 = 0; i1 < 3; i1++) {
        for (int i2 = (i1 + 1); i2 < 3; i2++) {
            auto mp2rage_filter = itk::BinaryFunctorImageFilter<QI::ImageXF, QI::ImageXF, QI::ImageXF, MP2Contrast<float>>::New();
            mp2rage_filter->SetInput1(vols[i1]->GetOutput());
            mp2rage_filter->SetInput2(vols[i2]->GetOutput());

            auto mask_filter = itk::MaskImageFilter<QI::ImageXF, QI::ImageF, QI::ImageXF>::New();
            mask_filter->SetInput1(mp2rage_filter->GetOutput());
            mask_filter->SetMaskImage(mask);
            mask_filter->Update();
            string outName = outPrefix + "MP3_C" + to_string(i1) + to_string(i2) + QI::OutExt();
            auto realFilter = itk::ComplexToRealImageFilter<QI::ImageXF, QI::ImageF>::New();
            realFilter->SetInput(mask_filter->GetOutput());
            realFilter->Update();
            conImages[ind] = realFilter->GetOutput();
            conImages[ind]->DisconnectPipeline();
            if (complex_output) {
                QI::writeResult<QI::ImageXF>(mask_filter->GetOutput(), outName);
            } else {
                QI::writeResult<QI::ImageF>(realFilter->GetOutput(), outName);
            }
            ind++;
        }
    }

    cout << "Attempting quantitative bit" << endl;
    auto apply = itk::ApplyAlgorithmFilter<MP3LookupAlgo>::New();
    auto lookup = make_shared<MP3LookupAlgo>();
    apply->SetAlgorithm(lookup);

    typedef itk::ComposeImageFilter<QI::ImageF, QI::VectorImageF> ComposeType;
    auto composer = ComposeType::New();
    composer->SetInput(0, conImages[0]);
    composer->SetInput(1, conImages[1]);
    composer->SetInput(2, conImages[2]);
    apply->SetInput(0, composer->GetOutput());
    apply->SetMask(mask);
    apply->Update();
    QI::writeResult(apply->GetOutput(0), outPrefix + "MP3_T1" + QI::OutExt());
    QI::writeResult(apply->GetOutput(1), outPrefix + "MP3_B1" + QI::OutExt());
    QI::writeResult(apply->GetOutput(2), outPrefix + "MP3_eta" + QI::OutExt());
    if (verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}

