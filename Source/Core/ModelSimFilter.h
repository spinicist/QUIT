/*
 *  ModelSimFilter.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_MODELSIMFILTER_H
#define QI_MODELSIMFILTER_H

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageToImageFilter.h"
#include "itkProgressReporter.h"
#include "itkTimeProbe.h"

#include "ImageTypes.h"
#include "ModelHelpers.h"

namespace QI {

template <typename ModelType, bool MultiOutput = false>
class ModelSimFilter
    : public itk::ImageToImageFilter<
          QI::VolumeF,
          itk::VectorImage<typename IOPrecision<typename ModelType::DataType>::Type, 3>> {
  public:
    static constexpr int ImageDim = 3;
    using OutputPixelType         = typename IOPrecision<typename ModelType::DataType>::Type;
    using OutputImageType         = itk::VectorImage<OutputPixelType, ImageDim>;
    using Self                    = ModelSimFilter;
    using Superclass              = itk::ImageToImageFilter<QI::VolumeF, OutputImageType>;
    using Pointer                 = itk::SmartPointer<Self>;
    using RegionType              = typename OutputImageType::RegionType;
    QI_ForwardNewMacro(Self);
    itkTypeMacro(Self, Superclass); /** Run-time type information (and related methods). */

    ModelSimFilter(ModelType const &m, bool const vb, std::string const &subregion) :
        m_model{m}, m_verbose{vb} {
        this->SetNumberOfRequiredInputs(ModelType::NV);
        if constexpr (MultiOutput) {
            this->SetNumberOfRequiredOutputs(m_model.num_outputs());
            for (size_t i = 0; i < m_model.num_outputs(); i++) {
                this->SetNthOutput(i, this->MakeOutput(i));
            }
        } else {
            this->SetNumberOfRequiredOutputs(1);
            this->SetNthOutput(0, this->MakeOutput(0));
        }
        if (subregion != "") {
            m_subregion    = RegionFromString<RegionType>(subregion);
            m_hasSubregion = true;
        }
    }

    void SetVarying(unsigned int i, const QI::VolumeF *image) {
        if (static_cast<int>(i) < ModelType::NV) {
            this->SetNthInput(i, const_cast<QI::VolumeF *>(image));
        } else {
            itkExceptionMacro("Requested varying input " << i << " does not exist ("
                                                         << ModelType::NV << " inputs)");
        }
    }

    typename QI::VolumeF::ConstPointer GetVarying(const int i) const {
        if (i < ModelType::NV) {
            return static_cast<const QI::VolumeF *>(this->itk::ProcessObject::GetInput(i));
        } else {
            itkExceptionMacro("Requested varying input " << i << " does not exist ("
                                                         << ModelType::NV << " inputs)");
        }
    }

    void SetFixed(const int i, const QI::VolumeF *image) {
        if (i < ModelType::NF) {
            this->SetNthInput(ModelType::NV + i, const_cast<QI::VolumeF *>(image));
        } else {
            itkExceptionMacro("Requested fixed input " << i << " does not exist (" << ModelType::NF
                                                       << " const inputs)");
        }
    }

    typename QI::VolumeF::ConstPointer GetFixed(const int i) const {
        if (i < ModelType::NF) {
            size_t index = ModelType::NV + i;
            return static_cast<const QI::VolumeF *>(this->itk::ProcessObject::GetInput(index));
        } else {
            itkExceptionMacro("Requested const input " << i << " does not exist (" << ModelType::NF
                                                       << " const inputs)");
        }
    }

    void SetMask(const QI::VolumeF *mask) {
        this->SetNthInput(ModelType::NV + ModelType::NF, const_cast<QI::VolumeF *>(mask));
    }

    typename QI::VolumeF::ConstPointer GetMask() const {
        return static_cast<const QI::VolumeF *>(
            this->itk::ProcessObject::GetInput(ModelType::NF + ModelType::NV));
    }

    OutputImageType *GetOutput(const int i) {
        return dynamic_cast<OutputImageType *>(this->itk::ProcessObject::GetOutput(i));
    }

    void SetSubregion(const RegionType &sr) {
        m_subregion    = sr;
        m_hasSubregion = true;
    }

    void SetNoise(const double s) { m_sigma = s; }

  private:
    ModelSimFilter(const Self &); // purposely not implemented
    void operator=(const Self &); // purposely not implemented

  protected:
    ModelType  m_model;
    double     m_sigma = 0.0;
    const bool m_verbose;
    bool       m_hasSubregion = false;
    RegionType m_subregion;

    ModelSimFilter() {}
    ~ModelSimFilter() {}

    itk::DataObject::Pointer
    MakeOutput(itk::ProcessObject::DataObjectPointerArraySizeType idx) override {
        itk::DataObject::Pointer output;
        if (idx < this->GetNumberOfRequiredOutputs()) {
            auto img = OutputImageType::New();
            output   = img;
        } else {
            QI::Fail("No output {}", idx);
        }
        return output.GetPointer();
    }

    void GenerateOutputInformation() override {
        Superclass::GenerateOutputInformation();
        const auto ip = this->GetInput(0);
        for (size_t i = 0; i < this->GetNumberOfRequiredOutputs(); i++) {
            const auto op = this->GetOutput(i);
            op->SetRegions(ip->GetLargestPossibleRegion());
            op->SetSpacing(ip->GetSpacing());
            op->SetOrigin(ip->GetOrigin());
            op->SetDirection(ip->GetDirection());
            if constexpr (MultiOutput) {
                op->SetNumberOfComponentsPerPixel(m_model.output_size(i));
            } else {
                op->SetNumberOfComponentsPerPixel(m_model.sequence.size());
            }
            op->Allocate(true);
        }
    }

    virtual void GenerateData() override {
        auto region = this->GetInput(0)->GetLargestPossibleRegion();
        if (m_hasSubregion) {
            if (region.IsInside(m_subregion)) {
                region = m_subregion;
            } else {
                itkExceptionMacro("Specified subregion is not entirely inside image.");
            }
        }

        Info(m_verbose, "Simulating...");
        this->GetMultiThreader()->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
        this->GetMultiThreader()->template ParallelizeImageRegion<ImageDim>(
            region,
            [this](const RegionType &outputRegion) {
                this->DynamicThreadedGenerateData(outputRegion);
            },
            this);
        Info(m_verbose, "Finished simulating.");
    }

    void DynamicThreadedGenerateData(const RegionType &region) override {
        std::vector<itk::ImageRegionConstIterator<QI::VolumeF>> varying_iters(ModelType::NV);
        for (int i = 0; i < ModelType::NV; i++) {
            varying_iters[i] =
                itk::ImageRegionConstIterator<QI::VolumeF>(this->GetInput(i), region);
        }
        std::vector<itk::ImageRegionConstIterator<QI::VolumeF>> fixed_iters(ModelType::NF);
        for (int i = 0; i < ModelType::NF; i++) {
            if (this->GetFixed(i)) {
                fixed_iters[i] =
                    itk::ImageRegionConstIterator<QI::VolumeF>(this->GetFixed(i), region);
            }
        }

        itk::ImageRegionConstIterator<QI::VolumeF> mask_iter;
        const auto                                 mask = this->GetMask();
        if (mask) {
            mask_iter = itk::ImageRegionConstIterator<QI::VolumeF>(mask, region);
        }

        std::vector<itk::ImageRegionIterator<OutputImageType>> output_iters(
            this->GetNumberOfRequiredOutputs());
        for (size_t i = 0; i < this->GetNumberOfRequiredOutputs(); i++) {
            output_iters[i] = itk::ImageRegionIterator<OutputImageType>(this->GetOutput(i), region);
        }

        while (!output_iters[0].IsAtEnd()) {
            if (!mask || mask_iter.Get()) {
                QI_ARRAYN(double, ModelType::NV) varying;
                for (int i = 0; i < ModelType::NV; i++) {
                    varying[i] = varying_iters[i].Get();
                }
                Eigen::ArrayXd fixed = m_model.fixed_defaults;
                for (int i = 0; i < ModelType::NF; i++) {
                    if (this->GetFixed(i)) {
                        fixed[i] = fixed_iters[i].Get();
                    }
                }

                if constexpr (MultiOutput) {
                    const auto signals = m_model.signals(varying, fixed);
                    for (size_t i = 0; i < signals.size(); i++) {
                        const auto output    = add_noise(signals[i], m_sigma);
                        const auto output_io = output.template cast<OutputPixelType>().eval();
                        itk::VariableLengthVector<OutputPixelType> data_out(output_io.data(),
                                                                            output_io.rows());
                        output_iters[i].Set(data_out);
                    }
                } else {
                    const auto signal    = m_model.signal(varying, fixed);
                    const auto output    = add_noise(signal, m_sigma);
                    const auto output_io = output.template cast<OutputPixelType>().eval();
                    itk::VariableLengthVector<OutputPixelType> data_out(output_io.data(),
                                                                        output_io.rows());
                    output_iters[0].Set(data_out);
                }
            }
            if (mask) {
                ++mask_iter;
            }
            for (int i = 0; i < ModelType::NV; i++) {
                ++varying_iters[i];
            }
            for (int i = 0; i < ModelType::NF; i++) {
                if (this->GetFixed(i)) {
                    ++fixed_iters[i];
                }
            }
            for (size_t i = 0; i < this->GetNumberOfRequiredOutputs(); i++) {
                ++output_iters[i];
            }
        }
    }
};

} // namespace QI

#endif