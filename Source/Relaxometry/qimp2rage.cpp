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

#include "itkBinaryFunctorImageFilter.h"
#include "itkExtractImageFilter.h"

#include "ApplyAlgorithmFilter.h"
#include "ImageTypes.h"
#include "Util.h"
#include "ImageIO.h"
#include "Args.h"
#include "MPRAGESequence.h"
#include "SequenceCereal.h"

template<class T> class MPRAGEFunctor {
public:
    MPRAGEFunctor() {}
    ~MPRAGEFunctor() {}
    bool operator!=(const MPRAGEFunctor &) const { return false; }
    bool operator==(const MPRAGEFunctor &other) const { return !(*this != other); }

    inline T operator()(const std::complex<T> &ti1, const std::complex<T> &ti2) const
    {
        const T a1 = abs(ti1);
        const T a2 = abs(ti2);
        return real((conj(ti1)*ti2)/(a1*a1 + a2*a2));
    }
};

namespace itk {

class MPRAGELookUpFilter : public ImageToImageFilter<QI::VolumeF, QI::VolumeF>
{
public:

protected:
    std::vector<float> m_T1, m_con;

public:
    /** Standard class typedefs. */
    typedef QI::VolumeF                        TImage;
    typedef MPRAGELookUpFilter                 Self;
    typedef ImageToImageFilter<TImage, TImage> Superclass;
    typedef SmartPointer<Self>                 Pointer;
    typedef typename TImage::RegionType        RegionType;

    itkNewMacro(Self); /** Method for creation through the object factory. */
    itkTypeMacro(MPRAGELookUpFilter, ImageToImageFilter); /** Run-time type information (and related methods). */

    void SetInput(const TImage *img) ITK_OVERRIDE {
        this->SetNthInput(0, const_cast<TImage*>(img));
    }

    void SetSequence(QI::MP2RAGESequence &sequence) {
        MPRAGEFunctor<double> con;
        m_T1.clear();
        m_con.clear();
        for (float T1 = 0.25; T1 < 4.3; T1 += 0.001) {
            Eigen::Array3d tp; tp << T1, 1.0, 1.0; // Fix B1 and eta
            Eigen::Array2cd sig = sequence.signal(1., T1, 1.0, 1.0);
            double c = con(sig[0], sig[1]);
            m_T1.push_back(T1);
            m_con.push_back(c);
            //cout << m_pars.back().transpose() << " : " << m_cons.back().transpose() << std::endl;
        }
        std::cout << "Lookup table has " << m_T1.size() << " entries" << std::endl;
    }

protected:
    MPRAGELookUpFilter() {
        this->SetNumberOfRequiredInputs(1);
        this->SetNumberOfRequiredOutputs(1);
        this->SetNthOutput(0, this->MakeOutput(0));
    }
    ~MPRAGELookUpFilter() {}

    DataObject::Pointer MakeOutput(ProcessObject::DataObjectPointerArraySizeType idx) override {
        //std::cout <<  __PRETTY_FUNCTION__ << std::endl;
        if (idx == 0) {
            DataObject::Pointer output = (TImage::New()).GetPointer();
            return output.GetPointer();
        } else {
            std::cerr << "No output " << idx << std::endl;
            return NULL;
        }
    }

    void ThreadedGenerateData(const RegionType &region, ThreadIdType threadId) ITK_OVERRIDE {
        //std::cout <<  __PRETTY_FUNCTION__ << std::endl;
        ImageRegionConstIterator<TImage> inputIter(this->GetInput(), region);
        ImageRegionIterator<TImage> outputIter(this->GetOutput(), region);

        while(!inputIter.IsAtEnd()) {
            const double ival = inputIter.Get();
            double best_distance = std::numeric_limits<double>::max();
            int best_index = 0;
            for (int i = m_T1.size(); i > 0; i--) {
                double distance = fabs(m_con[i] - ival);
                if (distance < best_distance) {
                    best_distance = distance;
                    best_index = i;
                }
            }
            //cout << "Best index " << best_index << " distance " << best_distance << " pars " << outputs.transpose() << " data " << data_inputs.transpose() << " cons" << m_cons[best_index].transpose() << std::endl;
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
    args::ArgumentParser parser("Calculates T1/B1 maps from MP2/3-RAGE data\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> input_path(parser, "INPUT FILE", "Path to complex MP-RAGE data");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "OUTPREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    QI::ParseArgs(parser, argc, argv);
    itk::MultiThreader::SetGlobalDefaultNumberOfThreads(threads.Get());

    if (verbose) std::cout << "Opening input file " << QI::CheckPos(input_path) << std::endl;
    auto inFile = QI::ReadImage<QI::SeriesXF>(QI::CheckPos(input_path));

    if (verbose) std::cout << "Combining MP2 contrasts" << std::endl;
    typedef itk::BinaryFunctorImageFilter<QI::VolumeXF, QI::VolumeXF, QI::VolumeF, MPRAGEFunctor<float>> MPRageContrastFilterType;
    int nti = inFile->GetLargestPossibleRegion().GetSize()[3];
    auto MPContrastFilter = MPRageContrastFilterType::New();
    typedef itk::ExtractImageFilter<QI::SeriesXF, QI::VolumeXF> ExtractType;
    auto vol_i = ExtractType::New();
    auto vol_j = ExtractType::New();
    vol_i->SetDirectionCollapseToSubmatrix();
    vol_j->SetDirectionCollapseToSubmatrix();
    QI::SeriesXF::RegionType region_i = inFile->GetLargestPossibleRegion();
    region_i.GetModifiableSize()[3] = 0;
    QI::SeriesXF::RegionType region_j = region_i;
    std::string outName = outarg ? outarg.Get() : QI::StripExt(input_path.Get());
    for (int i = 0; i < nti - 1; i++) {
        region_i.GetModifiableIndex()[3] = i;
        vol_i->SetInput(inFile);
        vol_i->SetExtractionRegion(region_i);
        vol_i->Update();
        MPContrastFilter->SetInput1(vol_i->GetOutput());
        for (int j = (i + 1); j < nti; j++) {
            region_j.GetModifiableIndex()[3] = j;
            vol_j->SetInput(inFile);
            vol_j->SetExtractionRegion(region_j);
            vol_j->Update();
            MPContrastFilter->SetInput2(vol_j->GetOutput());
            MPContrastFilter->Update();
            QI::WriteImage(MPContrastFilter->GetOutput(), outName + "_contrast" + QI::OutExt());
        }
    }
    std::cout << "Calculating T1" << std::endl;
    auto mp2rage_sequence = QI::ReadSequence<QI::MP2RAGESequence>(std::cin, verbose);
    auto apply = itk::MPRAGELookUpFilter::New();
    apply->SetSequence(mp2rage_sequence);
    apply->SetInput(MPContrastFilter->GetOutput());
    apply->Update();
    QI::WriteImage(apply->GetOutput(0), outName + "_R1" + QI::OutExt());
    if (verbose) std::cout << "Finished." << std::endl;
    return EXIT_SUCCESS;
}

