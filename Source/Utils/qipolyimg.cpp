/*
 *  qipolyimg.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include "Eigen/Dense"

#include "itkImageSource.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "ImageTypes.h"
#include "Util.h"
#include "Polynomial.h"
#include "Args.h"
#include "ImageIO.h"
#include "IO.h"
#include "EigenCereal.h"

namespace itk {

class PolynomialImage : public ImageSource<QI::VolumeF> {
public:
    typedef QI::VolumeF            TImage;
    typedef PolynomialImage        Self;
    typedef ImageSource<TImage>    Superclass;
    typedef SmartPointer<Self>     Pointer;

    itkNewMacro(Self);
    itkTypeMacro(Self, ImageSource);

    void SetReferenceImage(const SmartPointer<TImage> img) {
        m_reference = img;
    }

    void SetPolynomial(const QI::Polynomial<3> &p) { m_poly = p; }
    void SetMask(const TImage *mask) { this->SetNthInput(1, const_cast<TImage*>(mask)); }
    void SetCenter(const Eigen::Array3d &c) { m_center = c; }
    void SetScale(const double s) { m_scale = s; }
    typename TImage::ConstPointer GetMask() const { return static_cast<const TImage *>(this->ProcessObject::GetInput(1)); }
    
    void GenerateOutputInformation() ITK_OVERRIDE {
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
    Eigen::Array3d m_center = Eigen::Array3d::Zero();
    QI::Polynomial<3> m_poly;
    double m_scale = 1.0;

    PolynomialImage() {}
    ~PolynomialImage() {}
    void GenerateData() ITK_OVERRIDE {
        typename TImage::Pointer output = this->GetOutput();
        itk::ImageRegionIteratorWithIndex<TImage> imageIt(output, output->GetLargestPossibleRegion());
        imageIt.GoToBegin();
        const auto mask = this->GetMask();
        itk::ImageRegionConstIterator<TImage> maskIter;
        if (mask) {
            //if (m_verbose) std::cout << "Counting voxels in mask..." << std::endl;
            maskIter = itk::ImageRegionConstIterator<TImage>(mask, output->GetLargestPossibleRegion());
            maskIter.GoToBegin();
        }
        while(!imageIt.IsAtEnd()) {
            if (!mask || maskIter.Get()) {
                TImage::PointType p;
                m_reference->TransformIndexToPhysicalPoint(imageIt.GetIndex(), p);
                double val = m_poly.value((QI::Eigenify(p.GetVectorFromOrigin()) - m_center) / m_scale);
                imageIt.Set(val);
            } else {
                imageIt.Set(0);
            }
            ++imageIt;
            if (mask)
                ++maskIter;
        }
    }

private:
    PolynomialImage(const Self &);
    void operator=(const Self &);
};

} // End namespace itk

int main(int argc, char **argv) {
    Eigen::initParallel();
    args::ArgumentParser parser("Creates an image from polynomial coefficients, which are read from stdin.\n"
                                "\nhttp://github.com/spinicist/QUIT");
    args::Positional<std::string> ref_path(parser, "REFERENCE", "Reference image space to create the polynomial");
    args::Positional<std::string> out_path(parser, "OUTPUT", "Output image path");
    args::HelpFlag help(parser, "HELP", "Show this help message", {'h', "help"});
    args::Flag     verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::ValueFlag<int> threads(parser, "THREADS", "Use N threads (default=4, 0=hardware limit)", {'T', "threads"}, 4);
    args::ValueFlag<int> order(parser, "ORDER", "Specify the polynomial order (default 2)", {'o',"order"}, 2);
    args::ValueFlag<std::string> mask(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    QI::ParseArgs(parser, argc, argv, verbose);
    itk::MultiThreader::SetGlobalMaximumNumberOfThreads(threads.Get());
    if (verbose) std::cout << "Reading reference image " << QI::CheckPos(ref_path) << std::endl;
    QI::VolumeF::Pointer reference = QI::ReadImage(QI::CheckPos(ref_path));

    Eigen::Array3d center;
    double scale;
    Eigen::ArrayXd coeffs;
    {
        if (verbose) std::cout << "Reading polynomial" << std::endl;
        cereal::JSONInputArchive output(std::cin);
        QI::ReadCereal(output, "center", center);
        QI::ReadCereal(output, "scale", scale);
        QI::ReadCereal(output, "coeffs", coeffs);
    }
    if (verbose) std::cout << "Center point is: " << center.transpose() << "\nScale is: " << scale << std::endl;

    QI::Polynomial<3> poly(order.Get());
    if (coeffs.rows() != poly.nterms()) {
        QI_EXCEPTION("Require " + std::to_string(poly.nterms()) + " terms for " + std::to_string(order.Get()) + " order polynomial");
    }
    poly.setCoeffs(coeffs);
    if (verbose) std::cout << "Generating image" << std::endl;
    auto image = itk::PolynomialImage::New();
    image->SetReferenceImage(reference);
    image->SetPolynomial(poly);
    if (mask) image->SetMask(QI::ReadImage(mask.Get()));
    image->SetCenter(center);
    image->SetScale(scale);
    image->Update();
    QI::WriteImage(image->GetOutput(), out_path.Get());
    if (verbose) std::cout << "Finished." << std::endl;
    return EXIT_SUCCESS;
}
