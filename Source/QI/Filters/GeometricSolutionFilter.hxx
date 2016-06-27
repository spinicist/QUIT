/*
 *  GeometricSolutionFilter.h
 *
 *  Copyright (c) 2016 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_GEOMETRICSOLUTIONFILTER_HXX
#define QI_GEOMETRICSOLUTIONFILTER_HXX

namespace QI {

void GeometricSolutionFilter::SetRegularise(const RegEnum &r) { m_Regularise = r;}
const GeometricSolutionFilter::RegEnum &GeometricSolutionFilter::GetRegularise()       { return m_Regularise; }
void GeometricSolutionFilter::SetPhases(const size_t p) {
    if (p < 4)
        QI_EXCEPTION("Must have a minimum of 4 phase-cycling patterns.");
    if ((p % 2) != 0)
        QI_EXCEPTION("Number of phases must be even.");

    m_phases = p;
    m_lines = m_phases / 2;
    m_crossings = Choose(m_lines, 2);
    this->Modified();
}
void GeometricSolutionFilter::SetMask(const TMask *mask) { this->SetNthInput(1, const_cast<TMask*>(mask)); }
typename GeometricSolutionFilter::TMask::ConstPointer GeometricSolutionFilter::GetMask() const { return static_cast<const TMask *>(this->ProcessObject::GetInput(1)); }

void GeometricSolutionFilter::GenerateOutputInformation() ITK_OVERRIDE {
    Superclass::GenerateOutputInformation();
    if ((this->GetInput()->GetNumberOfComponentsPerPixel() % m_phases) != 0) {
        QI_EXCEPTION("Input size and number of phases do not match");
    }
    m_flips = (this->GetInput()->GetNumberOfComponentsPerPixel() / m_phases);
    auto op = this->GetOutput();
    op->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    op->SetNumberOfComponentsPerPixel(m_flips);
    op->Allocate();
}

GeometricSolutionFilter::GeometricSolutionFilter() {
    this->SetNumberOfRequiredInputs(1);
    this->SetNumberOfRequiredOutputs(1);
    this->SetNthOutput(0, this->MakeOutput(0));
    this->SetPhases(4);
}

void GeometricSolutionFilter::ThreadedGenerateData(const TIn::RegionType &region, itk::ThreadIdType threadId) ITK_OVERRIDE {
    //std::cout <<  __PRETTY_FUNCTION__ << endl;
    itk::ImageRegionConstIterator<TIn> inputIter(this->GetInput(), region);
    auto m = this->GetMask();
    itk::ImageRegionConstIterator<TMask> maskIter;
    if (m) {
        maskIter = itk::ImageRegionConstIterator<TMask>(m, region);
    }
    itk::ImageRegionIterator<TOut> outputIter(this->GetOutput(), region);

    while(!inputIter.IsAtEnd()) {
        if (!m || maskIter.Get()) {
            itk::VariableLengthVector<std::complex<float>> inputVector = inputIter.Get();
            Eigen::Map<const Eigen::ArrayXcf> allInput(inputVector.GetDataPointer(), m_phases);
            const Eigen::ArrayXcd a = allInput.head(m_lines).cast<std::complex<double>>();
            const Eigen::ArrayXcd b = allInput.tail(m_lines).cast<std::complex<double>>();

            std::complex<double> sum(0., 0.);
            for (size_t i = 0; i < m_lines; i++) {
                for (size_t j = i + 1; j < m_lines; j++) {
                    const std::complex<double> di = b[i] -  a[i], dj = b[j] - a[j];
                    const std::complex<double> ni(di.imag(), -di.real()), nj(dj.imag(), -dj.real());

                    const double mu = cdot(a[j] - a[i], nj) / cdot(di, nj);
                    const double nu = cdot(a[i] - a[j], ni) / cdot(dj, ni);
                    const double xi = 1.0 - pow(cdot(di, dj) / (abs(di)*abs(dj)), 2.0);

                    const std::complex<double> cs = (a[i] + a[j] + b[i] + b[j]) / 4.0;
                    const std::complex<double> gs = a[i] + mu * di;

                    switch (m_Regularise) {
                    case RegEnum::None: sum += gs; break;
                    case RegEnum::Magnitude:
                        if (norm(gs) < std::max({std::norm(a[i]), std::norm(a[j]), std::norm(b[i]), std::norm(b[j])})) {
                            sum += gs;
                        } else {
                            sum += cs;
                        }
                        break;
                    case RegEnum::Line:
                        if ((mu > -xi) && (mu < 1 + xi) && (nu > -xi) && (nu < 1 + xi)) {
                            sum += gs;
                        } else {
                            sum += cs;
                        }
                        break;
                    }
                }
            }
            outputIter.Set(static_cast<std::complex<float>>(sum / static_cast<double>(m_crossings)));
        }
        ++inputIter;
        if (m)
            ++maskIter;
        ++outputIter;
    }
}

}

#endif // QI_GEOMETRICSOLUTIONFILTER_HXX