/*
 *  qipolyfit.cpp
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, you can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include "Eigen/Dense"

#include "itkImageToImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageMomentsCalculator.h"
#include "Types.h"
#include "Util.h"
#include "ImageIO.h"
#include "Polynomial.h"
#include "Args.h"
#include "Fit.h"

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

    itkSetMacro(Robust, bool);
    itkGetMacro(Robust, bool);

    //void SetInput(const TImage *img) ITK_OVERRIDE      { this->SetNthInput(0, const_cast<TImage*>(img)); }
    //typename TImage::ConstPointer GetInput() const { return static_cast<const TImage *>(this->ProcessObject::GetInput(0)); }

    const QI::Polynomial<3> &GetPolynomial() const { return m_poly; } 
    void SetPolynomial(const QI::Polynomial<3> &p) { m_poly = p; }
    void SetMask(const TImage *mask) { this->SetNthInput(1, const_cast<TImage*>(mask)); }
    void SetCenter(const itk::Point<double, 3>& v) { m_center = v; }
    typename TImage::ConstPointer GetMask() const { return static_cast<const TImage *>(this->ProcessObject::GetInput(1)); }
    void GenerateOutputInformation() ITK_OVERRIDE {
        Superclass::GenerateOutputInformation();
        auto op = this->GetOutput();
        op->SetRegions(this->GetInput()->GetLargestPossibleRegion());
        op->Allocate();
    }

protected:
    QI::Polynomial<3> m_poly;
    itk::Point<double, 3> m_center;
    bool m_Robust;

    PolynomialFitImageFilter() {
        this->SetNumberOfRequiredInputs(1);
        m_center.Fill(0.0);
    }
    ~PolynomialFitImageFilter() {}

    void GenerateData() ITK_OVERRIDE {
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
        Eigen::VectorXd y(N);
        itk::ImageRegionConstIteratorWithIndex<TImage> imageIt(input,region);
        imageIt.GoToBegin();
        int yi = 0;
        while(!imageIt.IsAtEnd()) {
            if (!mask || maskIter.Get()) {
                TImage::PointType p, p2;
                input->TransformIndexToPhysicalPoint(imageIt.GetIndex(), p);
                p2 = p - m_center;
                Eigen::Vector3d ep(p2[0], p2[1], p2[2]);
                X.row(yi) = m_poly.terms(ep);
                y[yi] = imageIt.Get();
                ++yi;
            }
            ++imageIt;
            if (mask)
                ++maskIter;
        }
        Eigen::VectorXd b = m_Robust ? QI::RobustLeastSquares(X, y) : QI::LeastSquares(X, y);
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
    args::Flag                    print_terms(parser, "TERMS", "Print out the polynomial terms", {"print-terms"});
    args::Flag                    robust(parser, "ROBUST", "Use a robust (Huber) fit", {'r', "robust"});
    args::ValueFlag<int>          order(parser, "ORDER", "Specify the polynomial order (default 4)", {'o',"order"}, 4);
    args::ValueFlag<std::string>  mask_path(parser, "MASK", "Only process voxels within the mask", {'m', "mask"});

    QI::ParseArgs(parser, argc, argv);
    if (verbose) std::cout << "Reading input from: " << QI::CheckPos(input_path) << std::endl;
    auto input = QI::ReadImage(QI::CheckPos(input_path));
    auto fit = itk::PolynomialFitImageFilter::New();
    QI::Polynomial<3> poly(order.Get());
    fit->SetInput(input);
    fit->SetPolynomial(poly);
    fit->SetRobust(robust);
    itk::Point<double, 3> center; center.Fill(0.0);
    if (mask_path) {
        if (verbose) std::cout << "Reading mask from: " << mask_path.Get() << std::endl;
        auto mask_image = QI::ReadImage(mask_path.Get());
        fit->SetMask(mask_image);
        auto moments = itk::ImageMomentsCalculator<QI::VolumeF>::New();
        moments->SetImage(mask_image);
        moments->Compute();
        // ITK seems to put a minus sign on CoG
        if (verbose) std::cout << "Mask CoG is: " << -moments->GetCenterOfGravity() << std::endl;
        center = -moments->GetCenterOfGravity();
    }
    fit->SetCenter(center);
    fit->Update();
    // itk::Point does not have symmetric operator<< and operator>>
    std::cout << center[0] << " " << center[1] << " " << center[2] << std::endl;
    std::cout << fit->GetPolynomial().coeffs().transpose() << std::endl;
    if (print_terms)
        fit->GetPolynomial().print_terms();
    if (verbose) std::cout << "Finished." << std::endl;
    return EXIT_SUCCESS;
}
