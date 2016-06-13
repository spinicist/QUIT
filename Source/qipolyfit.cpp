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

#include <getopt.h>
#include <iostream>
#include "Eigen/Dense"

#include "itkImageToImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkMaskImageFilter.h"
#include "QI/Types.h"
#include "QI/Util.h"
#include "QI/Polynomial.h"

using namespace std;
using namespace Eigen;

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
        VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
        m_poly.setCoeffs(b);
    }

private:
    PolynomialFitImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
};

} // End namespace itk

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: qipolyfit [options] input \n\
\n\
Fits a 3D polynomial to a volume and prints the co-efficients to stdout\n\
\n\
Options:\n\
    --help, -h        : Print this message.\n\
    --verbose, -v     : Print more information.\n\
    --mask, -m file   : Mask input with specified file.\n\
    --order, -o N     : Specify the polynomial order (default 2)\n\
    --threads, -T N   : Use N threads (default=hardware limit).\n"
};

const struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"verbose", no_argument, 0, 'v'},
    {"mask", required_argument, 0, 'm'},
    {"order", required_argument, 0, 'o'},
    {"threads", required_argument, 0, 'T'},
    {0, 0, 0, 0}
};
const char *short_options = "hvm:o:T:";

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
    Eigen::initParallel();

    bool verbose = false;
    QI::VolumeF::Pointer mask = ITK_NULLPTR;
    int indexptr = 0, c, order = 2;
    while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
        switch (c) {
            case 'v': verbose = true; break;
            case 'm':
                if (verbose) cout << "Reading mask file " << optarg << endl;
                mask = QI::ReadImage(optarg);
                break;
            case 'o':
                order = stoi(optarg);
                if (verbose) cout << "Polynomical order is: " << order << endl;
                break;
            case 'T': itk::MultiThreader::SetGlobalMaximumNumberOfThreads(atoi(optarg)); break;
            case 'h':
                cout << QI::GetVersion() << endl << usage << endl;
                return EXIT_SUCCESS;
            case '?': // getopt will print an error message
                return EXIT_FAILURE;
            default:
                cout << "Unhandled option " << string(1, c) << endl;
                return EXIT_FAILURE;
        }
    }
    if ((argc - optind) != 1) {
        cout << "Incorrect number of arguments." << endl << usage << endl;
        return EXIT_FAILURE;
    }
    QI::VolumeF::Pointer input = QI::ReadImage(argv[optind]);
    auto fit = itk::PolynomialFitImageFilter::New();
    QI::Polynomial poly(order);
    fit->SetInput(input);
    fit->SetPolynomial(poly);
    fit->SetMask(mask);
    fit->Update();
    cout << fit->GetPolynomial().coeffs().transpose() << endl;
    if (verbose) cout << "Finished." << endl;
    return EXIT_SUCCESS;
}
