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

#include <Eigen/Core>

#include "Args.h"
#include "Fit.h"
#include "ImageIO.h"
#include "ImageTypes.h"
#include "JSON.h"
#include "Polynomial.h"
#include "Util.h"
#include "itkImageMomentsCalculator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageToImageFilter.h"

namespace itk {

class PolynomialFitImageFilter : public ImageToImageFilter<QI::VolumeF, QI::VolumeF> {
  public:
    /** Standard class typedefs. */
    typedef QI::VolumeF TImage;

    typedef PolynomialFitImageFilter           Self;
    typedef ImageToImageFilter<TImage, TImage> Superclass;
    typedef SmartPointer<Self>                 Pointer;
    typedef typename TImage::RegionType        RegionType;

    itkNewMacro(Self);
    itkTypeMacro(Self, Superclass);

    itkSetMacro(Robust, bool);
    itkGetMacro(Robust, bool);

    // void SetInput(const TImage *img) ITK_OVERRIDE      { this->SetNthInput(0,
    // const_cast<TImage*>(img)); } typename TImage::ConstPointer GetInput() const { return
    // static_cast<const TImage *>(this->ProcessObject::GetInput(0)); }

    const QI::Polynomial<3> &GetPolynomial() const { return m_poly; }
    void                     SetPolynomial(const QI::Polynomial<3> &p) { m_poly = p; }
    void SetMask(const TImage *mask) { this->SetNthInput(1, const_cast<TImage *>(mask)); }
    void SetCenter(const Eigen::Array3d &v) { m_center = v; }
    void SetScale(const double s) { m_scale = s; }
    typename TImage::ConstPointer GetMask() const {
        return static_cast<const TImage *>(this->ProcessObject::GetInput(1));
    }
    double GetResidual() { return m_residual; }
    void   GenerateOutputInformation() ITK_OVERRIDE { Superclass::GenerateOutputInformation(); }

  protected:
    QI::Polynomial<3> m_poly;
    Eigen::Array3d    m_center = Eigen::Array3d::Zero();
    bool              m_Robust = false;
    double            m_scale = 1.0, m_residual = 0.0;

    PolynomialFitImageFilter() {
        this->SetNumberOfRequiredInputs(1);
        this->SetNumberOfRequiredOutputs(0);
    }
    ~PolynomialFitImageFilter() {}

    void GenerateData() ITK_OVERRIDE {
        typename TImage::ConstPointer input  = this->GetInput();
        auto                          region = input->GetLargestPossibleRegion();

        const auto                       mask = this->GetMask();
        ImageRegionConstIterator<TImage> maskIter;
        int                              N = 0;
        if (mask) {
            // if (m_verbose) std::cout << "Counting voxels in mask..." );
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
        Eigen::MatrixXd                                X(N, m_poly.nterms());
        Eigen::VectorXd                                y(N);
        itk::ImageRegionConstIteratorWithIndex<TImage> imageIt(input, region);
        imageIt.GoToBegin();
        int yi = 0;
        while (!imageIt.IsAtEnd()) {
            if (!mask || maskIter.Get()) {
                TImage::PointType p;
                input->TransformIndexToPhysicalPoint(imageIt.GetIndex(), p);
                X.row(yi) =
                    m_poly.terms((QI::Eigenify(p.GetVectorFromOrigin()) - m_center) / m_scale);
                y[yi] = imageIt.Get();
                ++yi;
            }
            ++imageIt;
            if (mask)
                ++maskIter;
        }
        Eigen::VectorXd b = m_Robust ? QI::RobustLeastSquares(X, y, &m_residual) :
                                       QI::LeastSquares(X, y, &m_residual);
        m_poly.setCoeffs(b);
    }

  private:
    PolynomialFitImageFilter(const Self &); // purposely not implemented
    void operator=(const Self &);           // purposely not implemented
};

} // End namespace itk

int polyfit_main(args::Subparser &parser) {
    args::Positional<std::string> input_path(parser, "INPUT", "Input file");
    args::Flag print_terms(parser, "TERMS", "Print out the polynomial terms", {"print-terms"});
    args::Flag robust(parser, "ROBUST", "Use a robust (Huber) fit", {'r', "robust"});
    args::ValueFlag<int> order(
        parser, "ORDER", "Specify the polynomial order (default 4)", {'o', "order"}, 4);
    args::ValueFlag<std::string> mask_path(
        parser, "MASK", "Only process voxels within the mask", {'m', "mask"});
    parser.Parse();

    auto input = QI::ReadImage(QI::CheckPos(input_path), verbose);
    auto fit   = itk::PolynomialFitImageFilter::New();

    QI::Polynomial<3> poly(order.Get());
    fit->SetInput(input);
    fit->SetPolynomial(poly);
    fit->SetRobust(robust);
    Eigen::Array3d center = Eigen::Array3d::Zero();
    if (mask_path) {
        auto mask_image = QI::ReadImage(mask_path.Get(), verbose);
        fit->SetMask(mask_image);
        auto moments = itk::ImageMomentsCalculator<QI::VolumeF>::New();
        moments->SetImage(mask_image);
        moments->Compute();
        // ITK seems to put a minus sign on CoG
        QI::Log(verbose, "Mask CoG is: {}", -moments->GetCenterOfGravity());
        center = QI::Eigenify(-moments->GetCenterOfGravity());
    }
    fit->SetCenter(center);
    QI::VolumeF::SizeType    size    = input->GetLargestPossibleRegion().GetSize();
    QI::VolumeF::SpacingType spacing = input->GetSpacing();
    for (int i = 0; i < 3; i++)
        spacing[i] *= size[i];
    double scale = (spacing / 2).GetNorm();
    fit->SetScale(scale);
    fit->Update();
    json doc;
    doc["center"] = center;
    doc["scale"]  = scale;
    doc["coeffs"] = fit->GetPolynomial().coeffs();
    if (print_terms)
        doc["terms"] = fit->GetPolynomial().get_terms();
    QI::WriteJSON(std::cout, doc);
    QI::Log(verbose, "Residual: {}", fit->GetResidual());
    QI::Log(verbose, "Finished.");
    return EXIT_SUCCESS;
}
