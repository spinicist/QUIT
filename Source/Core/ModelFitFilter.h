/*
 *  ModelFitFilter.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_MODELFITFILTER_H
#define QI_MODELFITFILTER_H

#include <vector>
#include <tuple>
#include <Eigen/Core>

#include "itkImageToImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVariableLengthVector.h"
#include "itkVectorImage.h"
#include "itkTimeProbe.h"

#include "Macro.h"
#include "ModelHelpers.h"
#include "FitFunction.h"

namespace itk {

template<typename FitType>
class ModelFitFilter : public ImageToImageFilter<VectorImage<typename IOPrecision<typename FitType::InputType>::Type, 3>,
                                                 typename BlockTypes<FitType::Blocked, 3, typename IOPrecision<typename FitType::OutputType>::Type>::Type> {
public:
    static const int ImageDimension = 3;
    static const bool Blocked = FitType::Blocked;
    static const bool Indexed = FitType::Indexed;
    using InputType  = typename FitType::InputType;
    using OutputType = typename FitType::OutputType;
    using ParameterType  = typename FitType::ModelType::ParameterType;

    using InputPixelType  = typename IOPrecision<InputType>::Type;
    using OutputPixelType = typename IOPrecision<OutputType>::Type;
    using FixedPixelType  = typename IOPrecision<ParameterType>::Type;

    using TInputImage  = VectorImage<InputPixelType, ImageDimension>;
    using TFixedImage  = Image<FixedPixelType, ImageDimension>;
    using TMaskImage   = Image<float, ImageDimension>;

    using TOutputImage   = typename BlockTypes<Blocked, ImageDimension, OutputPixelType>::Type;
    using TFlagImage     = typename BlockTypes<Blocked, ImageDimension, typename FitType::FlagType>::Type;
    using TResidualImage = typename BlockTypes<Blocked, ImageDimension, InputPixelType>::Type;
    using TResidualsImage = VectorImage<InputPixelType, ImageDimension>;

    using TRegion    = typename TInputImage::RegionType;
    using TIndex     = typename TRegion::IndexType;

    using Self       = ModelFitFilter;
    using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
    using Pointer    = SmartPointer<Self>;

    QI_ForwardNewMacro(Self);
    itkTypeMacro(ModelFitFilter, ImageToImageFilter); /** Run-time type information (and related methods). */

    static const int FlagOutputOffset = 0;
    static const int ResidualOutputOffset = 1;
    static const int ResidualsOutputOffset = 2;

    ModelFitFilter(const FitType *f) :
        m_fit(f)
    {
        this->SetNumberOfRequiredInputs(m_fit->n_inputs());
        this->SetNumberOfRequiredOutputs(m_fit->n_outputs()+ResidualsOutputOffset);
        for (int i = 0; i < (m_fit->n_outputs()+ResidualsOutputOffset+m_fit->n_inputs()); i++) {
            this->SetNthOutput(i, this->MakeOutput(i));
        }
    }

    void SetInput(unsigned int i, const TInputImage *image) override {
        if (static_cast<int>(i) < m_fit->n_inputs()) {
            this->SetNthInput(i, const_cast<TInputImage*>(image));
        } else {
            itkExceptionMacro("Requested input " << i << " does not exist (" << m_fit->n_inputs() << " inputs)");
        }
    }

    typename TInputImage::ConstPointer GetInput(const int i) const {
        if (i < m_fit->n_inputs()) {
            return static_cast<const TInputImage *> (this->ProcessObject::GetInput(i));
        } else {
            QI_FAIL("Requested input " << i << " does not exist (" << m_fit->n_inputs() << " inputs)");
        }
    }

    void SetFixed(const int i, const TFixedImage *image) {
        if (i < m_fit->n_fixed()) {
            this->SetNthInput(m_fit->n_inputs() + 1 + i, const_cast<TFixedImage*>(image));
        } else {
            itkExceptionMacro("Requested const input " << i << " does not exist (" << m_fit->n_fixed() << " const inputs)");
        }
    }

    typename TFixedImage::ConstPointer GetFixed(const int i) const {
        if (i < m_fit->n_fixed()) {
            size_t index = m_fit->n_inputs() + 1 + i;
            return static_cast<const TFixedImage *> (this->ProcessObject::GetInput(index));
        } else {
            QI_FAIL("Requested const input " << i << " does not exist (" << m_fit->n_fixed() << " const inputs)");
        }
    }

    void SetMask(const TMaskImage *mask) {
        this->SetNthInput(m_fit->n_inputs(), const_cast<TMaskImage *>(mask));
    }

    typename TMaskImage::ConstPointer GetMask() const {
        return static_cast<const TMaskImage *>(this->ProcessObject::GetInput(m_fit->n_inputs()));
    }

    void SetOutputAllResiduals(const bool r) { m_allResiduals = r; }
    void SetVerbose(const bool v) { m_verbose = v; }
    
    void SetSubregion(const TRegion &sr) {
        if (m_verbose) std::cout << "Setting subregion to: " << std::endl << sr << std::endl;
        m_subregion = sr;
        m_hasSubregion = true;
    }

    void SetBlocks(const int &nb) {
        if constexpr(Blocked) {
            m_blocks = nb;
        } else {
            QI_FAIL("Cannot set block size on a non-blocked filter");
        }
    }

    TOutputImage *GetOutput(const int i) {
        if (i < m_fit->n_outputs()) {
            return dynamic_cast<TOutputImage *>(this->ProcessObject::GetOutput(i));
        } else {
            QI_FAIL("Requested output " << std::to_string(i) << " is past maximum (" << std::to_string(m_fit->n_outputs()) << ")");
        }
    }

    TResidualImage *GetResidualOutput() {
        return dynamic_cast<TResidualImage *>(this->ProcessObject::GetOutput(m_fit->n_outputs()+ResidualOutputOffset));
    }

    TResidualsImage *GetResidualsOutput(const int i) {
        return dynamic_cast<TResidualsImage *>(this->ProcessObject::GetOutput(m_fit->n_outputs()+ResidualsOutputOffset+i));
    }

    TFlagImage *GetFlagOutput() {
        return dynamic_cast<TFlagImage *>(this->ProcessObject::GetOutput(m_fit->n_outputs()+FlagOutputOffset));
    }

    RealTimeClock::TimeStampType GetTotalTime() const { return m_elapsedTime; }

private:
    ModelFitFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented

protected:
    DataObject::Pointer MakeOutput(ProcessObject::DataObjectPointerArraySizeType idx) override {
        DataObject::Pointer output;
        using itype = ProcessObject::DataObjectPointerArraySizeType; // Stop unsigned long versus int warnings
        if (idx < static_cast<itype>(m_fit->n_outputs())) {
            output = (TOutputImage::New()).GetPointer();
        } else if (idx == static_cast<itype>(m_fit->n_outputs() + FlagOutputOffset)) {
            auto img = TFlagImage::New();
            output = img;
        } else if (idx == static_cast<itype>(m_fit->n_outputs() + ResidualOutputOffset)) {
            auto img = TResidualImage::New();
            output = img;
        } else if (idx < static_cast<itype>(m_fit->n_outputs() + ResidualsOutputOffset + m_fit->n_inputs())) {
            auto img = TResidualImage::New();
            output = img;
        } else {
            QI_FAIL("Attempted to create output " << idx << ", index too high");
        }
        return output.GetPointer();
    }

    const FitType *m_fit;
    bool m_verbose = false, m_hasSubregion = false, m_allResiduals = false;
    TRegion m_subregion;
    int m_blocks = 1;

    RealTimeClock::TimeStampType m_elapsedTime = 0.0;

    virtual void GenerateOutputInformation() override {
        Superclass::GenerateOutputInformation();
        for (int i = 0; i < m_fit->n_inputs(); i++) {
            if ((m_fit->input_size(i) * m_blocks) != static_cast<int>(this->GetInput(i)->GetNumberOfComponentsPerPixel())) {
                QI_FAIL("Input " << i << " has incorrect size " << this->GetInput(i)->GetNumberOfComponentsPerPixel() << ", should be " << (m_fit->input_size(i) * m_blocks));
            }
        }

        auto input     = this->GetInput(0);
        auto region    = input->GetLargestPossibleRegion();
        auto spacing   = input->GetSpacing();
        auto origin    = input->GetOrigin();
        auto direction = input->GetDirection();
        if (m_verbose) std::cout << "Allocating output memory" << std::endl;
        for (int i = 0; i < m_fit->n_outputs(); i++) {
            auto op = this->GetOutput(i);
            op->SetRegions(region);
            op->SetSpacing(spacing);
            op->SetOrigin(origin);
            op->SetDirection(direction);
            if constexpr(Blocked) {
                op->SetNumberOfComponentsPerPixel(m_blocks);
            }
            op->Allocate(true);
        }
        if (m_allResiduals) {
            for (int i = 0; i < m_fit->n_inputs(); i++) {
                if (m_verbose) std::cout << "Allocating memory for point residuals " << i << std::endl;
                auto r = this->GetResidualsOutput(i);
                r->SetRegions(region);
                r->SetSpacing(spacing);
                r->SetOrigin(origin);
                r->SetDirection(direction);
                r->SetNumberOfComponentsPerPixel(m_fit->input_size(i) * m_blocks);
                r->Allocate(true);
            }
        }
        if (m_verbose) std::cout << "Allocating total residual memory" << std::endl;
        auto r = this->GetResidualOutput();
        r->SetRegions(region);
        r->SetSpacing(spacing);
        r->SetOrigin(origin);
        r->SetDirection(direction);
        if constexpr(Blocked) {
            r->SetNumberOfComponentsPerPixel(m_blocks);
        }
        r->Allocate(true);
        if (m_verbose) std::cout << "Allocating flag memory" << std::endl;
        auto f = this->GetFlagOutput();
        f->SetRegions(region);
        f->SetSpacing(spacing);
        f->SetOrigin(origin);
        f->SetDirection(direction);
        if constexpr(Blocked) {
            f->SetNumberOfComponentsPerPixel(m_blocks);
        }
        f->Allocate(true);
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

        this->GetMultiThreader()->template ParallelizeImageRegion<ImageDimension>(
            region,
            [this](const typename TOutputImage::RegionType &outputRegion) { this->DynamicThreadedGenerateData(outputRegion); },
            this);
        // DynamicThreadedGenerateData(region);
    }

    virtual void DynamicThreadedGenerateData(const TRegion &region) override {
        ImageRegionConstIterator<TMaskImage> mask_iter;
        const auto mask = this->GetMask();
        if (mask) {
            mask_iter = ImageRegionConstIterator<TMaskImage>(mask, region);
        }
        
        std::vector<ImageRegionConstIterator<TInputImage>> input_iters(m_fit->n_inputs());
        std::vector<ImageRegionIterator<TResidualsImage>> residuals_iters(m_fit->n_inputs());
        for (int i = 0; i < m_fit->n_inputs(); i++) {
            input_iters[i] = ImageRegionConstIterator<TInputImage>(this->GetInput(i), region);
            if (m_allResiduals) {
                residuals_iters[i] = ImageRegionIterator<TResidualsImage>(this->GetResidualsOutput(i), region);
            }
        }
        std::vector<ImageRegionConstIterator<TFixedImage>> fixed_iters(m_fit->n_fixed());
        for (int i = 0; i < m_fit->n_fixed(); i++) {
            typename TFixedImage::ConstPointer c = this->GetFixed(i);
            if (c) {
                fixed_iters[i] = ImageRegionConstIterator<TFixedImage>(c, region);
            }
        }

        std::vector<ImageRegionIterator<TOutputImage>> output_iters(m_fit->n_outputs());
        for (int i = 0; i < m_fit->n_outputs(); i++) {
            output_iters[i] = ImageRegionIterator<TOutputImage>(this->GetOutput(i), region);
        }

        // Keep the index on this one to report voxel locations where algorithm fails.
        ImageRegionIteratorWithIndex<TResidualImage> residual_iter(this->GetResidualOutput(), region);
        ImageRegionIterator<TFlagImage> flag_iter(this->GetFlagOutput(), region);

        using InputArray    = QI_ARRAY(typename FitType::InputType);
        using OutputArray   = QI_ARRAYN(typename FitType::OutputType, FitType::ModelType::NV);
        using ResidualArray = QI_ARRAY(typename FitType::ResidualType);
        while(!input_iters[0].IsAtEnd()) {
            if (!mask || mask_iter.Get()) {
                for (int b = 0; b < m_blocks; b++) {
                    std::vector<InputArray> inputs(m_fit->n_inputs());
                    for (int i = 0; i < m_fit->n_inputs(); i++) {
                        auto input_data = input_iters[i].Get();
                        inputs[i] = InputArray(m_fit->input_size(i));
                        const int block_start = b * m_fit->input_size(i);
                        for (Eigen::Index j = 0; j < m_fit->input_size(i); j++) {
                            inputs[i][j] = input_data[j + block_start];
                        }
                    }

                    auto fixed = m_fit->model.fixed_defaults;
                    for (size_t i = 0; i < fixed_iters.size(); i++) {
                        if (this->GetFixed(i)) {
                            fixed[i] = fixed_iters[i].Get();
                        }
                    }

                    OutputArray outputs;
                    typename FitType::ResidualType residual{};
                    typename FitType::FlagType flag{};
                    std::vector<ResidualArray> residuals; // Leave size 0 if user doesn't want them
                    if (m_allResiduals) {
                        for (int i = 0; i < m_fit->n_inputs(); i++) {
                            const int block_size = m_fit->input_size(i);
                            residuals.push_back(ResidualArray::Zero(block_size));
                        }
                    }

                    QI::FitReturnType status;
                    if constexpr(Blocked && Indexed) {
                        status = m_fit->fit(inputs, fixed, outputs, residual, residuals, flag, b, residual_iter.GetIndex());
                    } else if constexpr(Blocked) { 
                        status = m_fit->fit(inputs, fixed, outputs, residual, residuals, flag, b);
                    } else if constexpr(Indexed) {
                        status = m_fit->fit(inputs, fixed, outputs, residual, residuals, flag, residual_iter.GetIndex());
                    } else {
                        status = m_fit->fit(inputs, fixed, outputs, residual, residuals, flag);
                    }

                    if (!std::get<0>(status)) {
                        std::cerr << "Fit failed for voxel " << residual_iter.GetIndex() << ": " << std::get<1>(status) << std::endl;
                    }

                    if constexpr(Blocked) {
                        flag_iter.Get()[b] = flag;
                        residual_iter.Get()[b] = residual;
                        for (int i = 0; i < m_fit->n_outputs(); i++) {
                            output_iters[i].Get()[b] = outputs[i];
                        }
                    } else {
                        flag_iter.Set(flag);
                        residual_iter.Set(static_cast<InputType>(residual));
                        for (int i = 0; i < m_fit->n_outputs(); i++) {
                        output_iters[i].Set(outputs[i]);
                        }
                    }
                    if (m_allResiduals) {
                        for (int i = 0; i < m_fit->n_inputs(); i++) {
                            const int block_start = m_fit->input_size(i) * b;
                            for (int j = 0; j < m_fit->input_size(i); j++) {
                                residuals_iters[i].Get()[j + block_start] = residuals[i][j];
                            }
                        }
                    }
                }
            } else {
                if constexpr(Blocked) {
                    flag_iter.Get().Fill(0);
                    residual_iter.Get().Fill(0.0f);
                    for (int i = 0; i < m_fit->n_outputs(); i++) {
                    output_iters[i].Get().Fill(0.0f);
                    }
                } else {
                    flag_iter.Set(0);
                    flag_iter.Set(0.0f);
                    for (int i = 0; i < m_fit->n_outputs(); i++) {
                    output_iters[i].Set(0.0f);
                    }
                }
                if (m_allResiduals) {
                    for (int i = 0; i < m_fit->n_inputs(); i++) {
                        residuals_iters[i].Get().Fill(0.0f);
                    }
                }
            }
            
            if (this->GetMask())
                ++mask_iter;
            for (int i = 0; i < m_fit->n_inputs(); i++) {
                ++input_iters[i];
                if (m_allResiduals)
                    ++residuals_iters[i];
            }
            for (int i = 0; i < m_fit->n_fixed(); i++) {
                if (this->GetFixed(i))
                    ++fixed_iters[i];
            }
            for (int i = 0; i < m_fit->n_outputs(); i++) {
                ++output_iters[i];
            }
            ++flag_iter;
            ++residual_iter;
        }
    }

};

} // End namespace QI

#endif // MODELFITFILTER_H
