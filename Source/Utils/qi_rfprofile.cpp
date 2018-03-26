/*
 *  qi_rfprofile.cpp
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include "Eigen/Dense"
#include <unsupported/Eigen/Splines>
#include <cereal/types/vector.hpp>

#include "itkImageSource.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkMaskImageFilter.h"
#include "ImageTypes.h"
#include "Util.h"
#include "Args.h"
#include "ImageIO.h"
#include "IO.h"
#include "Spline.h"
#include "EigenCereal.h"

namespace itk {

class ProfileImage : public ImageSource<QI::VectorVolumeF> {
public:
    using Self       = ProfileImage;
    using Superclass = ImageSource<QI::VectorVolumeF>;
    using TRegion    = QI::VectorVolumeF::RegionType;
    using Pointer    = SmartPointer<Self>;

    itkNewMacro(Self);
    itkTypeMacro(Self, ImageSource);

    void SetReference(const SmartPointer<QI::VolumeF> ref) {
        m_reference = ref;
    }

    void SetRF(const Eigen::ArrayXd pos, const std::vector<Eigen::ArrayXd> vals) {
        m_splines.clear();
        for (const auto &v : vals) {
            m_splines.push_back(QI::SplineInterpolator(pos, v));
        }
    }

    void SetMask(const QI::VolumeF *mask) { this->SetNthInput(1, const_cast<QI::VolumeF*>(mask)); }
    typename QI::VolumeF::ConstPointer GetMask() const { return static_cast<const QI::VolumeF *>(this->ProcessObject::GetInput(1)); }
    
    void GenerateOutputInformation() ITK_OVERRIDE {
        auto output = this->GetOutput();
        output->SetNumberOfComponentsPerPixel(m_splines.size());
        output->SetRegions(m_reference->GetLargestPossibleRegion());
        output->SetSpacing(m_reference->GetSpacing());
        output->SetDirection(m_reference->GetDirection());
        output->SetOrigin(m_reference->GetOrigin());
        output->Allocate();
    }

protected:
    SmartPointer<QI::VolumeF> m_reference;
    std::vector<QI::SplineInterpolator> m_splines;

    ProfileImage(){
    }
    ~ProfileImage(){}
    void ThreadedGenerateData(const TRegion &region, ThreadIdType threadId) ITK_OVERRIDE {
        auto output = this->GetOutput();
        ImageSliceIteratorWithIndex<QI::VectorVolumeF> imageIt(output, region);
        imageIt.SetFirstDirection(0);
        imageIt.SetSecondDirection(1);
        imageIt.GoToBegin();
        ImageSliceIteratorWithIndex<QI::VolumeF> it_B1(m_reference, region);
        it_B1.SetFirstDirection(0);
        it_B1.SetSecondDirection(1);
        it_B1.GoToBegin();

        const auto mask = this->GetMask();
        ImageSliceConstIteratorWithIndex<QI::VolumeF> maskIter;
        if (mask) {
            maskIter = ImageSliceConstIteratorWithIndex<QI::VolumeF>(mask, region);
            maskIter.SetFirstDirection(0);
            maskIter.SetSecondDirection(1);
            maskIter.GoToBegin();
        }

        // Calculate geometric center
        QI::VectorVolumeF::IndexType idx_center;
        for (int i = 0; i < 3; i++) {
            idx_center[i] = m_reference->GetLargestPossibleRegion().GetSize()[i] / 2;
        }
        QI::VectorVolumeF::PointType pt_center; m_reference->TransformIndexToPhysicalPoint(idx_center, pt_center);
        itk::VariableLengthVector<float> zero(m_splines.size()); zero.Fill(0.);

        ProgressReporter progress(this, threadId, region.GetNumberOfPixels(), 10);
        while(!imageIt.IsAtEnd()) {
            QI::VectorVolumeF::PointType pt, pt_rf;
            m_reference->TransformIndexToPhysicalPoint(imageIt.GetIndex(), pt);
            pt_rf = pt - pt_center;
            itk::VariableLengthVector<float> vals(m_splines.size());

            for (int i = 0; i < m_splines.size(); i++) {
                vals[i] = m_splines[i](pt_rf[2]);
            }
            while (!imageIt.IsAtEndOfSlice()) {
                while (!imageIt.IsAtEndOfLine()) {
                    if (!mask || maskIter.Get()) {
                        imageIt.Set(vals * it_B1.Get());
                    } else {
                        imageIt.Set(zero);
                    }
                    ++imageIt;
                    ++it_B1;
                    if (mask) {
                        ++maskIter;
                    }
                    if (threadId == 0) {
                        progress.CompletedPixel();
                    }
                }
                imageIt.NextLine();
                it_B1.NextLine();
                if (mask)
                    maskIter.NextLine();
            }
            imageIt.NextSlice();
            it_B1.NextSlice();
            if (mask) {
                maskIter.NextSlice();
            }
        }
    }

private:
    ProfileImage(const Self &);
    void operator=(const Self &);
};

} // End namespace itk

int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Generates a relative B1 map from a B1+ and RF profile.\nInput is the B1+ map.\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> b1plus_path(parser, "B1+_FILE", "Input B1+ file");
    args::Positional<std::string> output_path(parser, "B1_FILE", "Output relative B1 file");
    
    args::HelpFlag help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag     debug(parser, "DEBUG", "Output debugging messages", {'d', "debug"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "PREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> subregion(parser, "REGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    QI::ParseArgs(parser, argc, argv, verbose);

    itk::MultiThreader::SetGlobalMaximumNumberOfThreads(threads.Get());
    if (verbose) std::cout << "Reading image " << QI::CheckPos(b1plus_path) << std::endl;
    auto reference = QI::ReadImage(QI::CheckPos(b1plus_path));

    cereal::JSONInputArchive input(std::cin);
    Eigen::ArrayXd rf_pos;
    std::vector<Eigen::ArrayXd> rf_vals;
    if (verbose) std::cout << "Reading slab profile" << std::endl;
    QI::ReadCereal(input, "rf_pos", rf_pos);
    QI::ReadCereal(input, "rf_vals", rf_vals);
    try {
        input(cereal::make_nvp("rf_vals", rf_vals));
    } catch (cereal::RapidJSONException &e) {
        std::cerr << "Could not read array parameter 'rf_vals'" << std::endl;
        return EXIT_FAILURE;
    }
    for (const auto &v : rf_vals) {
        if (rf_pos.rows() != v.rows()) {
            QI_FAIL("Number of points must match number of values");
        };
    }
    if (verbose) {
        std::cout << "Profile has " << rf_pos.rows() << " points." << std::endl;
        if (verbose) std::cout << "Generating image" << std::endl;
    }
    auto image = itk::ProfileImage::New();
    image->SetReference(reference);
    image->SetRF(rf_pos, rf_vals);
    if (mask) image->SetMask(QI::ReadImage(mask.Get()));
    auto monitor = QI::GenericMonitor::New();
    image->AddObserver(itk::ProgressEvent(), monitor);
    image->Update();
    if (verbose) std::cout << "Finished, writing output: " << QI::CheckPos(output_path) << std::endl;
    QI::WriteVectorImage(image->GetOutput(), QI::CheckPos(output_path));
    return EXIT_SUCCESS;
}
