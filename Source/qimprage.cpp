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
"A tool to process MP-2/3-RAGE images\n\
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

template<class T> class MPRAGEFunctor {
public:
    MPRAGEFunctor() {}
    ~MPRAGEFunctor() {}
    bool operator!=(const MPRAGEFunctor &) const { return false; }
    bool operator==(const MPRAGEFunctor &other) const { return !(*this != other); }

    inline T operator()(const complex<T> &ti1, const complex<T> &ti2) const
    {
        const T a1 = abs(ti1);
        const T a2 = abs(ti2);
        return real((conj(ti1)*ti2)/(a1*a1 + a2*a2));
    }
};

namespace itk {

class MPRAGELookUpFilter : public ImageToImageFilter<QI::ImageF, QI::ImageF>
{
public:

protected:
    std::vector<float> m_T1, m_con;

public:
    /** Standard class typedefs. */
    typedef QI::ImageF                         TImage;
    typedef MPRAGELookUpFilter                 Self;
    typedef ImageToImageFilter<TImage, TImage> Superclass;
    typedef SmartPointer<Self>                 Pointer;
    typedef typename TImage::RegionType        RegionType;

    itkNewMacro(Self); /** Method for creation through the object factory. */
    itkTypeMacro(MPRAGELookUpFilter, ImageToImageFilter); /** Run-time type information (and related methods). */

    void SetInput(const TImage *img) override {
        this->SetNthInput(0, const_cast<TImage*>(img));
    }

protected:
    MPRAGELookUpFilter() {
        this->SetNumberOfRequiredInputs(1);
        this->SetNumberOfRequiredOutputs(1);
        this->SetNthOutput(0, this->MakeOutput(0));
        MP2RAGE sequence(true);
        MPRAGEFunctor<double> con;
        m_T1.clear();
        m_con.clear();
        for (float T1 = 0.5; T1 < 4.3; T1 += 0.001) {
            Array3d tp; tp << T1, 1.0, 1.0; // Fix B1 and eta
            Array2cd sig = sequence.signal(1., T1, 1.0, 1.0);
            double c = con(sig[0], sig[1]);
            m_T1.push_back(T1);
            m_con.push_back(c);
            //cout << m_pars.back().transpose() << " : " << m_cons.back().transpose() << endl;
        }
        cout << "Lookup table has " << m_T1.size() << " entries" << endl;
    }
    ~MPRAGELookUpFilter() {}

    DataObject::Pointer MakeOutput(unsigned int idx) {
        //std::cout <<  __PRETTY_FUNCTION__ << endl;
        if (idx == 0) {
            DataObject::Pointer output = (TImage::New()).GetPointer();
            return output.GetPointer();
        } else {
            std::cerr << "No output " << idx << std::endl;
            return NULL;
        }
    }

    virtual void ThreadedGenerateData(const RegionType &region, ThreadIdType threadId) override {
        //std::cout <<  __PRETTY_FUNCTION__ << endl;
        ImageRegionConstIterator<TImage> inputIter(this->GetInput(), region);
        ImageRegionIterator<TImage> outputIter(this->GetOutput(), region);

        while(!inputIter.IsAtEnd()) {
            const double ival = inputIter.Get();
            double best_distance = numeric_limits<double>::max();
            int best_index = 0;
            for (int i = 0; i < m_T1.size(); i++) {
                double distance = fabs(m_con[i] - ival);
                if (distance < best_distance) {
                    best_distance = distance;
                    best_index = i;
                }
            }
            //cout << "Best index " << best_index << " distance " << best_distance << " pars " << outputs.transpose() << " data " << data_inputs.transpose() << " cons" << m_cons[best_index].transpose() << endl;
            outputIter.Set(1./m_T1[best_index]);
            ++inputIter;
            ++outputIter;
        }
    }

private:
    MPRAGELookUpFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
};

} // End namespace itk

int main(int argc, char **argv) {
    int indexptr = 0, c;
    string outName = "";
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
                outName = optarg;
                cout << "Output prefix will be: " << outName << endl;
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
    string inName(argv[optind]);
    if (verbose) cout << "Opening input file " << inName << endl;
    if (outName == "")
        outName = QI::StripExt(inName) + "_";
    auto inFile = QI::ReadTimeseriesXF::New();
    inFile->SetFileName(inName);
    inFile->Update();
    if (verbose) cout << "Processing" << endl;

    typedef itk::BinaryFunctorImageFilter<QI::ImageXF, QI::ImageXF, QI::ImageF, MPRAGEFunctor<float>> MPRageContrastFilterType;
    int nti = inFile->GetOutput()->GetLargestPossibleRegion().GetSize()[3];

    auto MPContrastFilter = MPRageContrastFilterType::New();
    typedef itk::ExtractImageFilter<QI::TimeseriesXF, QI::ImageXF> ExtractType;
    auto vol_i = ExtractType::New();
    auto vol_j = ExtractType::New();
    vol_i->SetDirectionCollapseToSubmatrix();
    vol_j->SetDirectionCollapseToSubmatrix();
    QI::TimeseriesXF::RegionType region_i = inFile->GetOutput()->GetLargestPossibleRegion();
    region_i.GetModifiableSize()[3] = 0;
    QI::TimeseriesXF::RegionType region_j = region_i;
    for (int i = 0; i < nti - 1; i++) {
        region_i.GetModifiableIndex()[3] = i;
        vol_i->SetInput(inFile->GetOutput());
        vol_i->SetExtractionRegion(region_i);
        vol_i->Update();
        MPContrastFilter->SetInput1(vol_i->GetOutput());
        for (int j = (i + 1); j < nti; j++) {
            region_j.GetModifiableIndex()[3] = j;
            vol_j->SetInput(inFile->GetOutput());
            vol_j->SetExtractionRegion(region_j);
            vol_j->Update();
            MPContrastFilter->SetInput2(vol_j->GetOutput());
            MPContrastFilter->Update();
            QI::writeResult(MPContrastFilter->GetOutput(), outName + "_contrast" + QI::OutExt());
        }
    }

    if (do_lookup) {
        auto apply = itk::MPRAGELookUpFilter::New();
        apply->SetInput(MPContrastFilter->GetOutput());
        apply->Update();
        QI::writeResult(apply->GetOutput(0), outName + "_R1" + QI::OutExt());
    }
    if (verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}

