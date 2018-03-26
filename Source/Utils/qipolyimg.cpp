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
#include "itkMaskImageFilter.h"
#include "ImageTypes.h"
#include "Util.h"
#include "Polynomial.h"
#include "Args.h"
#include "ImageIO.h"
#include "IO.h"

using namespace std;
using namespace Eigen;

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
    void SetCenter(const itk::Point<double, 3>& v) { m_center = v; }
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
    itk::Point<double, 3> m_center;
    QI::Polynomial<3> m_poly;

    PolynomialImage(){
        m_center.Fill(0.0);
    }
    ~PolynomialImage(){}
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
                TImage::PointType p, p2;
                m_reference->TransformIndexToPhysicalPoint(imageIt.GetIndex(), p);
                p2 = p - m_center;
                Eigen::Vector3d ep(p2[0], p2[1], p2[2]);
                double val = m_poly.value(ep);
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
    if (verbose) cout << "Reading reference image " << QI::CheckPos(ref_path) << std::endl;
    QI::VolumeF::Pointer reference = QI::ReadImage(QI::CheckPos(ref_path));
    itk::Point<double, 3> center;
    std::cin >> center;
    if (verbose) std::cout << "Center point is: " << center << std::endl;
    if (verbose) cout << "Building polynomial" << std::endl;
    QI::Polynomial<3> poly(order.Get());
    ArrayXd coeff;
    std::string dummy; std::getline(cin, dummy); // Damn C++ stream operators
    QI::ReadArray(cin, coeff);
    if (coeff.rows() != poly.nterms()) {
        QI_EXCEPTION("Require " + to_string(poly.nterms()) + " terms for " + to_string(order.Get()) + " order polynomial");
    }
    poly.setCoeffs(coeff);
    if (verbose) cout << "Generating image" << std::endl;
    auto image = itk::PolynomialImage::New();
    image->SetReferenceImage(reference);
    image->SetPolynomial(poly);
    if (mask) image->SetMask(QI::ReadImage(mask.Get()));
    image->SetCenter(center);
    image->Update();
    QI::WriteImage(image->GetOutput(), out_path.Get());
    if (verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}
