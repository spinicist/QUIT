/*
 *  qipolyfit.cpp
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

#include "itkImageToImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkMaskImageFilter.h"
#include "QI/Types.h"
#include "QI/Util.h"
#include "QI/IO.h"
#include "QI/Polynomial.h"
#include "QI/Args.h"

namespace itk {

class PolynomialFitImageFilter : public ImageToImageFilter<QI::VolumeF, QI::VolumeF> {
public:
    /** Standard class typedefs. */
    typedef QI::VolumeF     TImage;

    typedef PolynomialFitImageFilter           Self;
    typedef ImageToImageFilter<TImage, TImage> Superclass;
    typedef SmartPointer<Self>                 Pointer;
    typedef typename TImage::RegionType        RegionType;

    itkNewMacro(Self);
    itkTypeMacro(Self, Superclass);

    //void SetInput(const TImage *img) ITK_OVERRIDE      { this->SetNthInput(0, const_cast<TImage*>(img)); }
    //typename TImage::ConstPointer GetInput() const { return static_cast<const TImage *>(this->ProcessObject::GetInput(0)); }

    const QI::Polynomial &GetPolynomial() const { return m_poly; } 
    void SetPolynomial(const QI::Polynomial &p) { m_poly = p; }
    void SetMask(const TImage *mask) { this->SetNthInput(1, const_cast<TImage*>(mask)); }
    typename TImage::ConstPointer GetMask() const { return static_cast<const TImage *>(this->ProcessObject::GetInput(1)); }
    virtual void GenerateOutputInformation() ITK_OVERRIDE {
        Superclass::GenerateOutputInformation();
        auto op = this->GetOutput();
        op->SetRegions(this->GetInput()->GetLargestPossibleRegion());
        op->Allocate();
    }

protected:
    QI::Polynomial m_poly;
    PolynomialFitImageFilter() {
        this->SetNumberOfRequiredInputs(1);
    }
    ~PolynomialFitImageFilter() {}

    virtual void GenerateData() ITK_OVERRIDE {
        typename TImage::ConstPointer input = this->GetInput();
        auto region = input->GetLargestPossibleRegion();

        const auto mask = this->GetMask();
        ImageRegionConstIterator<TImage> maskIter;
        int N = 0;
        if (mask) {
            //if (m_verbose) std::cout << "Counting voxels in mask..." << std::endl;
            maskIter = ImageRegionConstIterator<TImage>(mask, region);
            maskIter.GoToBegin();
            while (!maskIter.IsAtEnd()) {
                if (maskIter.Get())
                    ++N;
                ++maskIter;
            }
            maskIter.GoToBegin(); // Reset
        } else {
            N = region.GetNumberOfPixels();
        }
        Eigen::MatrixXd X(N, m_poly.nterms());
        Eigen::VectorXd Y(N);
        itk::ImageRegionConstIteratorWithIndex<TImage> imageIt(input,region);
        imageIt.GoToBegin();
        ++imageIt;
        int yi = 0;
        while(!imageIt.IsAtEnd()) {
            if (!mask || maskIter.Get()) {
                TImage::PointType p;
                input->TransformIndexToPhysicalPoint(imageIt.GetIndex(), p);
                Eigen::Vector3d ep(p[0], p[1], p[2]);
                X.row(yi) = m_poly.terms(ep);
                Y[yi] = imageIt.Get();
                ++yi;
            }
            ++imageIt;
            if (mask)
                ++maskIter;
        }
        Eigen::VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
        m_poly.setCoeffs(b);
    }

private:
    PolynomialFitImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
};

} // End namespace itk

int main(int argc, char **argv) {
    Eigen::initParallel();

    args::ArgumentParser parser("Fits a 3D polynomial to a volume and prints the co-efficients to stdout.\n"
                                "http://github.com/spinicist/QUIT");

    args::Positional<std::string> input_path(parser, "INPUT", "Input file");
    args::HelpFlag                help(parser, "HELP", "Show this help menu", {'h', "help"});
    args::Flag                    verbose(parser, "VERBOSE", "Print more information", {'v', "verbose"});
    args::Flag                    print_terms(parser, "TERMS","Print out the polynomial terms", {"print-terms"});
    args::ValueFlag<int>          order(parser, "ORDER", "Specify the polynomial order (default 4)", {'o',"order"}, 4);
    args::ValueFlag<std::string>  mask_path(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});

    QI::ParseArgs(parser, argc, argv);
    if (verbose) std::cout << "Reading input from: " << QI::CheckPos(input_path) << std::endl;
    auto input = QI::ReadImage(QI::CheckPos(input_path));
    auto fit = itk::PolynomialFitImageFilter::New();
    QI::Polynomial poly(order.Get());
    fit->SetInput(input);
    fit->SetPolynomial(poly);
    if (mask_path) {
        if (verbose) std::cout << "Reading mask from: " << mask_path.Get() << std::endl;
        fit->SetMask(QI::ReadImage(mask_path.Get()));
    }
    fit->Update();
    std::cout << fit->GetPolynomial().coeffs().transpose() << std::endl;
    if (print_terms)
        fit->GetPolynomial().print_terms();
    if (verbose) std::cout << "Finished." << std::endl;
    return EXIT_SUCCESS;
}
