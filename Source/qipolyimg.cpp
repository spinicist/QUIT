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
#include "itkMaskImageFilter.h"
#include "QI/Types.h"
#include "QI/Util.h"
#include "QI/Polynomial.h"
#include "QI/Option.h"

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

    void SetPolynomial(const QI::Polynomial &p) { m_poly = p; }
    void SetMask(const TImage *mask) { this->SetNthInput(1, const_cast<TImage*>(mask)); }
    void SetCenter(const itk::Point<double, 3>& v) { m_center = v; }
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
    itk::Point<double, 3> m_center;
    QI::Polynomial m_poly;

    PolynomialImage(){
        m_center.Fill(0.0);
    }
    ~PolynomialImage(){}
    virtual void GenerateData() ITK_OVERRIDE {
        typename TImage::Pointer output = this->GetOutput();
        itk::ImageRegionIteratorWithIndex<TImage> imageIt(output, output->GetLargestPossibleRegion());
        imageIt.GoToBegin();
        const auto mask = this->GetMask();
        ImageRegionConstIterator<TImage> maskIter;
        if (mask) {
            //if (m_verbose) std::cout << "Counting voxels in mask..." << std::endl;
            maskIter = ImageRegionConstIterator<TImage>(mask, output->GetLargestPossibleRegion());
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
    QI::OptionList opts("Usage is: qipolyimg [options] reference output\n\nPolynomial coefficients are read from stdin.\n");
    QI::Option<int> num_threads(4,'T',"threads","Use N threads (default=4, 0=hardware limit)", opts);
    QI::Option<int> order(2,'o',"order","Specify the polynomial order (default 2)", opts);
    QI::ImageOption<QI::VolumeF> mask('m', "mask", "Mask input with specified file", opts);
    QI::Switch verbose('v',"verbose","Print more information", opts);
    QI::Help help(opts);
    std::deque<std::string> nonopts = opts.parse(argc, argv);
    if (nonopts.size() != 2) {
        std::cerr << opts << std::endl;
        std::cerr << "Required inputs are reference image and output filename." << std::endl;
        return EXIT_FAILURE;
    }
    itk::MultiThreader::SetGlobalMaximumNumberOfThreads(*num_threads);
    if (*verbose) cout << "Reading image " << argv[optind] << std::endl;
    QI::VolumeF::Pointer reference = QI::ReadImage(nonopts[0]);
    itk::Point<double, 3> center;
    std::cin >> center;
    if (*verbose) std::cout << "Center point is: " << center << std::endl;
    if (*verbose) cout << "Building polynomial" << std::endl;
    QI::Polynomial poly(*order);
    ArrayXd coeff;
    std::string dummy; std::getline(cin, dummy); // Damn C++ stream operators
    QI::ReadArray(cin, coeff);
    if (coeff.rows() != poly.nterms()) {
        QI_EXCEPTION("Require " + to_string(poly.nterms()) + " terms for " + to_string(*order) + " order polynomial");
    }
    poly.setCoeffs(coeff);
    if (*verbose) cout << "Generating image" << std::endl;
    auto image = itk::PolynomialImage::New();
    image->SetReferenceImage(reference);
    image->SetPolynomial(poly);
    image->SetMask(*mask);
    image->SetCenter(center);
    image->Update();
    QI::WriteImage(image->GetOutput(), nonopts[1]);
    if (*verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}
