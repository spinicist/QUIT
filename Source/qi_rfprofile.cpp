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

#include "itkImageSource.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkMaskImageFilter.h"
#include "QI/Types.h"
#include "QI/Util.h"
#include "QI/Args.h"
#include "QI/IO.h"
#include "QI/Fit.h"

namespace itk {

class ProfileImage : public ImageSource<QI::VolumeF> {
public:
    typedef QI::VolumeF         TImage;
    typedef ProfileImage        Self;
    typedef ImageSource<TImage> Superclass;
    typedef SmartPointer<Self>  Pointer;
    typedef typename TImage::RegionType TRegion;

    itkNewMacro(Self);
    itkTypeMacro(Self, ImageSource);

    void SetReferenceImage(const SmartPointer<TImage> img) {
        m_reference = img;
    }

    void SetProfile(const Eigen::ArrayXd &pos, const Eigen::ArrayXd &val) {
        m_spline = QI::SplineInterpolator(pos, val);
    }

    void SetMask(const TImage *mask) { this->SetNthInput(1, const_cast<TImage*>(mask)); }
    typename TImage::ConstPointer GetMask() const { return static_cast<const TImage *>(this->ProcessObject::GetInput(1)); }
    
    virtual void GenerateOutputInformation() ITK_OVERRIDE {
        Superclass::GenerateOutputInformation();
        auto output = this->GetOutput();
        output->SetRegions(m_reference->GetLargestPossibleRegion());
        output->SetSpacing(m_reference->GetSpacing());
        output->SetDirection(m_reference->GetDirection());
        output->SetOrigin(m_reference->GetOrigin());
        output->Allocate();
    }

protected:
    SmartPointer<TImage> m_reference;
    QI::SplineInterpolator m_spline;

    ProfileImage(){
    }
    ~ProfileImage(){}
    virtual void ThreadedGenerateData(const TRegion &region, ThreadIdType threadId) ITK_OVERRIDE {
        typename TImage::Pointer output = this->GetOutput();
        ImageSliceIteratorWithIndex<TImage> imageIt(output, region);
        imageIt.SetFirstDirection(0);
        imageIt.SetSecondDirection(1);
        imageIt.GoToBegin();
        ImageSliceIteratorWithIndex<TImage> it_B1(m_reference, region);
        it_B1.SetFirstDirection(0);
        it_B1.SetSecondDirection(1);
        it_B1.GoToBegin();

        const auto mask = this->GetMask();
        ImageSliceConstIteratorWithIndex<TImage> maskIter;
        if (mask) {
            maskIter = ImageSliceConstIteratorWithIndex<TImage>(mask, region);
            maskIter.SetFirstDirection(0);
            maskIter.SetSecondDirection(1);
            maskIter.GoToBegin();
        }

        // Calculate geometric center
        TImage::IndexType idx_center;
        for (int i = 0; i < 3; i++) {
            idx_center[i] = m_reference->GetLargestPossibleRegion().GetSize()[i] / 2;
        }
        TImage::PointType pt_center; m_reference->TransformIndexToPhysicalPoint(idx_center, pt_center);
        ProgressReporter progress(this, threadId, region.GetNumberOfPixels(), 10);
        while(!imageIt.IsAtEnd()) {
            TImage::PointType pt, pt_rf;
            m_reference->TransformIndexToPhysicalPoint(imageIt.GetIndex(), pt);
            pt_rf = pt - pt_center;
            const double val = m_spline(pt_rf[2]);
            while (!imageIt.IsAtEndOfSlice()) {
                while (!imageIt.IsAtEndOfLine()) {
                    if (!mask || maskIter.Get()) {
                        imageIt.Set(val * it_B1.Get());
                    } else {
                        imageIt.Set(0);
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
    args::Flag     noprompt(parser, "NOPROMPT", "Suppress input prompts", {'n', "no-prompt"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<std::string> outarg(parser, "PREFIX", "Add a prefix to output filenames", {'o', "out"});
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    args::ValueFlag<std::string> subregion(parser, "REGION", "Process subregion starting at voxel I,J,K with size SI,SJ,SK", {'s', "subregion"});
    QI::ParseArgs(parser, argc, argv);
    bool prompt = !noprompt;

    itk::MultiThreader::SetGlobalMaximumNumberOfThreads(threads.Get());
    if (verbose) std::cout << "Reading image " << QI::CheckPos(b1plus_path) << std::endl;
    auto reference = QI::ReadImage(QI::CheckPos(b1plus_path));

    if (prompt) std::cout << "Enter RF profile z-locations:" << std::endl;
    Eigen::ArrayXd rf_pos; QI::ReadArray(std::cin, rf_pos);
    if (prompt) std::cout << "Enter RF profile values:" << std::endl;
    Eigen::ArrayXd rf_val; QI::ReadArray(std::cin, rf_val);
    if (rf_pos.rows() != rf_val.rows()) {
        QI_FAIL("Number of points must match number of values");
    } else if (verbose) {
        std::cout << "Profile has " << rf_pos.rows() << " points." << std::endl;
    }
    
    if (verbose) std::cout << "Generating image" << std::endl;
    auto image = itk::ProfileImage::New();
    image->SetReferenceImage(reference);
    image->SetProfile(rf_pos, rf_val);
    if (mask) image->SetMask(QI::ReadImage(mask.Get()));
    auto monitor = QI::GenericMonitor::New();
    image->AddObserver(itk::ProgressEvent(), monitor);
    image->Update();
    if (verbose) std::cout << "Finished, writing output: " << QI::CheckPos(output_path) << std::endl;
    QI::WriteImage(image->GetOutput(), QI::CheckPos(output_path));
    return EXIT_SUCCESS;
}
