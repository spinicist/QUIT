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

#include <Eigen/Core>
#include <array>
#include <functional>
#include <tuple>
#include <vector>

#include "itkCommand.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageToImageFilter.h"
#include "itkTimeProbe.h"
#include "itkVariableLengthVector.h"
#include "itkVectorImage.h"

#include "FitFunction.h"
#include "Log.h"
#include "ModelHelpers.h"
#include "Monitor.h"
#include "Util.h"

namespace QI {

template <typename FitType>
class ModelFitFilter
    : public itk::ImageToImageFilter<
          itk::VectorImage<typename IOPrecision<typename FitType::InputType>::Type, 3>,
          typename BlockTypes<FitType::Blocked,
                              3,
                              typename IOPrecision<typename FitType::OutputType>::Type>::Type> {
  public:
    using ModelType     = typename FitType::ModelType;
    using InputType     = typename FitType::InputType;
    using OutputType    = typename FitType::OutputType;
    using ParameterType = typename ModelType::ParameterType;

    using InputPixelType  = typename IOPrecision<InputType>::Type;
    using OutputPixelType = typename IOPrecision<OutputType>::Type;
    using FixedPixelType  = typename IOPrecision<ParameterType>::Type;

    static constexpr int ImageDim = 3;

    using TInputImage = itk::VectorImage<InputPixelType, ImageDim>;
    using TFixedImage = itk::Image<FixedPixelType, ImageDim>;
    using TMaskImage  = itk::Image<float, ImageDim>;

    static constexpr bool Blocked = FitType::Blocked;

    using TOutputImage   = typename BlockTypes<Blocked, ImageDim, OutputPixelType>::Type;
    using TFlagImage     = typename BlockTypes<Blocked, ImageDim, typename FitType::FlagType>::Type;
    using TResidualImage = typename BlockTypes<Blocked, ImageDim, InputPixelType>::Type;
    using TResidualsImage = itk::VectorImage<InputPixelType, ImageDim>;

    using TRegion = typename TInputImage::RegionType;
    using TIndex  = typename TRegion::IndexType;

    using Self       = ModelFitFilter;
    using Superclass = itk::ImageToImageFilter<TInputImage, TOutputImage>;
    using Pointer    = itk::SmartPointer<Self>;

    static constexpr bool Indexed    = FitType::Indexed;
    static constexpr bool HasDerived = ModelType::ND > 0;

    QI_ForwardNewMacro(Self);
    itkTypeMacro(ModelFitFilter,
                 ImageToImageFilter); /** Run-time type information (and related methods). */

    // Input image offsets
    static constexpr int MaskOffset  = 0;
    static constexpr int FixedOffset = MaskOffset + 1;

    // Output image offsets
    static constexpr int DerivedOutputOffset = HasDerived ? ModelType::NV : -1;
    static constexpr int FlagOutputOffset =
        HasDerived ? DerivedOutputOffset + ModelType::ND : ModelType::NV;
    static constexpr int ResidualOutputOffset  = FlagOutputOffset + 1;
    static constexpr int ResidualsOutputOffset = ResidualOutputOffset + 1;
    static constexpr int TotalOutputs          = ResidualsOutputOffset + 1;

    ModelFitFilter(FitType const *    f,
                   const bool         verbose,
                   const bool         allResids,
                   std::string const &subregion) :
        m_fit(f),
        m_verbose(verbose), m_allResiduals(allResids) {
        this->SetNumberOfRequiredInputs(m_fit->n_inputs());
        this->SetNumberOfRequiredOutputs(TotalOutputs);
        for (int i = 0; i < TotalOutputs; i++) {
            this->SetNthOutput(i, this->MakeOutput(i));
        }
        if (m_verbose) {
            auto monitor = QI::GenericMonitor::New();
            this->AddObserver(itk::ProgressEvent(), monitor);
        }
        if (subregion != "") {
            m_subregion    = RegionFromString<TRegion>(subregion);
            m_hasSubregion = true;
        }
    }

    void SetInput(unsigned int i, const TInputImage *image) override {
        if (static_cast<int>(i) < m_fit->n_inputs()) {
            this->SetNthInput(i, const_cast<TInputImage *>(image));
        } else {
            itkExceptionMacro("Requested input " << i << " does not exist (" << m_fit->n_inputs()
                                                 << " inputs)");
        }
    }

    typename TInputImage::ConstPointer GetInput(const int i) const {
        if (i < m_fit->n_inputs()) {
            return static_cast<const TInputImage *>(this->itk::ProcessObject::GetInput(i));
        } else {
            QI::Fail(
                "Requested input {} but {} has {}", i, typeid(FitType).name(), m_fit->n_inputs());
        }
    }

    void SetFixed(const int i, const TFixedImage *image) {
        if (i < ModelType::NF) {
            this->SetNthInput(m_fit->n_inputs() + FixedOffset + i,
                              const_cast<TFixedImage *>(image));
        } else {
            QI::Fail("Tried to set fixed input {} but {} has {}",
                     i,
                     typeid(ModelType).name(),
                     ModelType::NF);
        }
    }

    typename TFixedImage::ConstPointer GetFixed(const int i) const {
        if (i < ModelType::NF) {
            size_t index = m_fit->n_inputs() + FixedOffset + i;
            return static_cast<const TFixedImage *>(this->itk::ProcessObject::GetInput(index));
        } else {
            QI::Fail("Requested fixed input {} but {} has {}",
                     i,
                     typeid(ModelType).name(),
                     ModelType::NF);
        }
    }

    void SetMask(const TMaskImage *mask) {
        this->SetNthInput(m_fit->n_inputs() + MaskOffset, const_cast<TMaskImage *>(mask));
    }

    typename TMaskImage::ConstPointer GetMask() const {
        return static_cast<const TMaskImage *>(
            this->itk::ProcessObject::GetInput(m_fit->n_inputs() + MaskOffset));
    }

    void SetOutputAllResiduals(const bool r) { m_allResiduals = r; }

    void SetSubregion(const TRegion &sr) {
        if (m_verbose)
            std::cout << "Setting subregion to: " << std::endl << sr << std::endl;
        m_subregion    = sr;
        m_hasSubregion = true;
    }

    void SetBlocks(const int &nb) {
        if constexpr (Blocked) {
            m_blocks = nb;
        } else {
            QI::Fail("Cannot set block size on a non-blocked filter");
        }
    }

    TOutputImage *GetOutput(const int i) {
        if (i < ModelType::NV) {
            return dynamic_cast<TOutputImage *>(this->itk::ProcessObject::GetOutput(i));
        } else {
            QI::Fail("Requested varying output {} but {} has {}",
                     i,
                     typeid(ModelType).name(),
                     ModelType::NV);
        }
    }

    TResidualImage *GetResidualOutput() {
        return dynamic_cast<TResidualImage *>(
            this->itk::ProcessObject::GetOutput(ResidualOutputOffset));
    }

    TResidualsImage *GetResidualsOutput(const int i) {
        return dynamic_cast<TResidualsImage *>(
            this->itk::ProcessObject::GetOutput(ResidualsOutputOffset + i));
    }

    TFlagImage *GetFlagOutput() {
        return dynamic_cast<TFlagImage *>(this->itk::ProcessObject::GetOutput(FlagOutputOffset));
    }

    TOutputImage *GetDerivedOutput(const int i) {
        if constexpr (HasDerived) {
            if (i < ModelType::ND) {
                return dynamic_cast<TOutputImage *>(
                    this->itk::ProcessObject::GetOutput(DerivedOutputOffset + i));
            } else {
                QI::Fail("Requested derived output {} but {} has {}",
                         i,
                         typeid(ModelType).name(),
                         ModelType::ND);
            }
        } else {
            QI::Fail("No derived outputs for this model");
        }
    }

    /*
     * I would prefer this to go in the constructor, but the ForwardNew macro can't cope with
     * container construction, i.e. {} in the call itself.
     */
    void ReadInputs(std::vector<std::string> const &              inputs,
                    std::array<std::string, ModelType::NF> const &fixed,
                    std::string const &                           mask) {
        if (static_cast<size_t>(m_fit->n_inputs()) != inputs.size()) {
            QI::Fail("Number of input file paths did not match number of inputs for model");
        }

        for (int i = 0; i < m_fit->n_inputs(); i++) {
            SetInput(i, QI::ReadImage<TInputImage>(inputs[i], m_verbose));
        }
        for (int f = 0; f < ModelType::NF; f++) {
            if (fixed[f] != "")
                SetFixed(f, QI::ReadImage<TFixedImage>(fixed[f], m_verbose));
        }
        if (mask != "")
            SetMask(QI::ReadImage<TMaskImage>(mask, m_verbose));
    }
    void WriteOutputs(std::string const &prefix) {
        for (int i = 0; i < ModelType::NV; i++) {
            QI::WriteImage(
                GetOutput(i), prefix + m_fit->model.varying_names.at(i) + QI::OutExt(), m_verbose);
        }
        if constexpr (ModelType::ND > 0) {
            for (int i = 0; i < ModelType::ND; i++) {
                QI::WriteImage(GetDerivedOutput(i),
                               prefix + m_fit->model.derived_names.at(i) + QI::OutExt(),
                               m_verbose);
            }
        }
        QI::WriteImage(GetResidualOutput(), prefix + "SoS_residual" + QI::OutExt(), m_verbose);
        QI::WriteImage(GetFlagOutput(), prefix + "iterations" + QI::OutExt(), m_verbose);
        if (m_allResiduals) {
            for (int i = 0; i < m_fit->n_inputs(); i++) {
                QI::WriteImage(GetResidualsOutput(0),
                               prefix + "residuals_" + std::to_string(i) + QI::OutExt(),
                               m_verbose);
            }
        }
    }

  private:
    ModelFitFilter(const Self &); // purposely not implemented
    void operator=(const Self &); // purposely not implemented

  protected:
    itk::DataObject::Pointer
    MakeOutput(itk::ProcessObject::DataObjectPointerArraySizeType idx) override {
        using itype = itk::ProcessObject::DataObjectPointerArraySizeType; // Stop unsigned long
                                                                          // versus int warnings
        if (idx < static_cast<itype>(ModelType::NV)) {
            return TOutputImage::New().GetPointer();
        } else if (idx < static_cast<itype>(FlagOutputOffset)) {
            return TOutputImage::New().GetPointer();
        } else if (idx == static_cast<itype>(FlagOutputOffset)) {
            return TFlagImage::New().GetPointer();
        } else if (idx == static_cast<itype>(ResidualOutputOffset)) {
            return TResidualImage::New().GetPointer();
        } else if (idx < static_cast<itype>(ResidualsOutputOffset + m_fit->n_inputs())) {
            return TResidualImage::New().GetPointer();
        } else {
            QI::Fail("Attempted to create output {} but {} has {}",
                     idx,
                     typeid(FitType).name(),
                     TotalOutputs);
        }
    }

    const FitType *m_fit;
    const bool     m_verbose, m_allResiduals;
    bool           m_hasSubregion = false;
    TRegion        m_subregion;
    int            m_blocks = 1;

    virtual void GenerateOutputInformation() override {
        Superclass::GenerateOutputInformation();
        for (int i = 0; i < m_fit->n_inputs(); i++) {
            if ((m_fit->input_size(i) * m_blocks) !=
                static_cast<int>(this->GetInput(i)->GetNumberOfComponentsPerPixel())) {
                QI::Fail("Input {} has incorrect number of volumes {}, should be {}",
                         i,
                         this->GetInput(i)->GetNumberOfComponentsPerPixel(),
                         (m_fit->input_size(i) * m_blocks));
            }
        }

        Log(m_verbose, "Allocating output image memory");
        auto input     = this->GetInput(0);
        auto region    = input->GetLargestPossibleRegion();
        auto spacing   = input->GetSpacing();
        auto origin    = input->GetOrigin();
        auto direction = input->GetDirection();
        for (int i = 0; i < ModelType::NV; i++) {
            auto op = this->GetOutput(i);
            op->SetRegions(region);
            op->SetSpacing(spacing);
            op->SetOrigin(origin);
            op->SetDirection(direction);
            if constexpr (Blocked) {
                op->SetNumberOfComponentsPerPixel(m_blocks);
            }
            op->Allocate(true);
        }
        if constexpr (HasDerived) {
            for (int i = 0; i < ModelType::ND; i++) {
                auto op = this->GetDerivedOutput(i);
                op->SetRegions(region);
                op->SetSpacing(spacing);
                op->SetOrigin(origin);
                op->SetDirection(direction);
                if constexpr (Blocked) {
                    op->SetNumberOfComponentsPerPixel(m_blocks);
                }
                op->Allocate(true);
            }
        }
        if (m_allResiduals) {
            for (int i = 0; i < m_fit->n_inputs(); i++) {
                auto r = this->GetResidualsOutput(i);
                r->SetRegions(region);
                r->SetSpacing(spacing);
                r->SetOrigin(origin);
                r->SetDirection(direction);
                r->SetNumberOfComponentsPerPixel(m_fit->input_size(i) * m_blocks);
                r->Allocate(true);
            }
        }
        auto r = this->GetResidualOutput();
        r->SetRegions(region);
        r->SetSpacing(spacing);
        r->SetOrigin(origin);
        r->SetDirection(direction);
        if constexpr (Blocked) {
            r->SetNumberOfComponentsPerPixel(m_blocks);
        }
        r->Allocate(true);
        auto f = this->GetFlagOutput();
        f->SetRegions(region);
        f->SetSpacing(spacing);
        f->SetOrigin(origin);
        f->SetDirection(direction);
        if constexpr (Blocked) {
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

        Info(m_verbose, "Processing...");
        this->GetMultiThreader()->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
        this->GetMultiThreader()->template ParallelizeImageRegion<ImageDim>(
            region,
            [this](const typename TOutputImage::RegionType &outputRegion) {
                this->DynamicThreadedGenerateData(outputRegion);
            },
            this);
        Info(m_verbose, "Finished processing.");
    }

    virtual void DynamicThreadedGenerateData(const TRegion &region) override {
        itk::ImageRegionConstIterator<TMaskImage> mask_iter;
        const auto                                mask = this->GetMask();
        if (mask) {
            mask_iter = itk::ImageRegionConstIterator<TMaskImage>(mask, region);
        }

        std::vector<itk::ImageRegionConstIterator<TInputImage>> input_iters(m_fit->n_inputs());
        std::vector<itk::ImageRegionIterator<TResidualsImage>>  residuals_iters(m_fit->n_inputs());
        for (int i = 0; i < m_fit->n_inputs(); i++) {
            input_iters[i] = itk::ImageRegionConstIterator<TInputImage>(this->GetInput(i), region);
            if (m_allResiduals) {
                residuals_iters[i] =
                    itk::ImageRegionIterator<TResidualsImage>(this->GetResidualsOutput(i), region);
            }
        }
        std::vector<itk::ImageRegionConstIterator<TFixedImage>> fixed_iters(ModelType::NF);
        for (int i = 0; i < ModelType::NF; i++) {
            typename TFixedImage::ConstPointer c = this->GetFixed(i);
            if (c) {
                fixed_iters[i] = itk::ImageRegionConstIterator<TFixedImage>(c, region);
            }
        }

        std::vector<itk::ImageRegionIterator<TOutputImage>> output_iters(ModelType::NV);
        for (int i = 0; i < ModelType::NV; i++) {
            output_iters[i] = itk::ImageRegionIterator<TOutputImage>(this->GetOutput(i), region);
        }

        std::vector<itk::ImageRegionIterator<TOutputImage>> derived_iters(ModelType::ND);
        if constexpr (HasDerived) {
            for (int i = 0; i < ModelType::ND; i++) {
                derived_iters[i] =
                    itk::ImageRegionIterator<TOutputImage>(this->GetDerivedOutput(i), region);
            }
        }

        // Keep the index on this one to report voxel locations where algorithm fails.
        itk::ImageRegionIteratorWithIndex<TResidualImage> residual_iter(this->GetResidualOutput(),
                                                                        region);
        itk::ImageRegionIterator<TFlagImage>              flag_iter(this->GetFlagOutput(), region);

        using InputArray    = QI_ARRAY(typename FitType::InputType);
        using OutputArray   = QI_ARRAYN(typename FitType::OutputType, ModelType::NV);
        using ResidualArray = QI_ARRAY(typename FitType::ResidualType);
        while (!input_iters[0].IsAtEnd()) {
            if (!mask || mask_iter.Get()) {
                for (int b = 0; b < m_blocks; b++) {
                    std::vector<InputArray> inputs(m_fit->n_inputs());
                    for (int i = 0; i < m_fit->n_inputs(); i++) {
                        auto input_data       = input_iters[i].Get();
                        inputs[i]             = InputArray(m_fit->input_size(i));
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

                    OutputArray                    outputs;
                    typename FitType::ResidualType residual{};
                    typename FitType::FlagType     flag{};
                    std::vector<ResidualArray> residuals; // Leave size 0 if user doesn't want them
                    if (m_allResiduals) {
                        for (int i = 0; i < m_fit->n_inputs(); i++) {
                            const int block_size = m_fit->input_size(i);
                            residuals.push_back(ResidualArray::Zero(block_size));
                        }
                    }

                    QI::FitReturnType status;
                    if constexpr (Blocked && Indexed) {
                        status = m_fit->fit(inputs,
                                            fixed,
                                            outputs,
                                            residual,
                                            residuals,
                                            flag,
                                            b,
                                            residual_iter.GetIndex());
                    } else if constexpr (Blocked) {
                        status = m_fit->fit(inputs, fixed, outputs, residual, residuals, flag, b);
                    } else if constexpr (Indexed) {
                        status = m_fit->fit(inputs,
                                            fixed,
                                            outputs,
                                            residual,
                                            residuals,
                                            flag,
                                            residual_iter.GetIndex());
                    } else {
                        status = m_fit->fit(inputs, fixed, outputs, residual, residuals, flag);
                    }

                    if (!status.success) {
                        QI::Warn("Fit failed for voxel {}: {}",
                                 residual_iter.GetIndex(),
                                 status.message);
                    }

                    if constexpr (Blocked) {
                        flag_iter.Get()[b]     = flag;
                        residual_iter.Get()[b] = residual;
                        for (int i = 0; i < ModelType::NV; i++) {
                            output_iters[i].Get()[b] = outputs[i];
                        }
                    } else {
                        flag_iter.Set(flag);
                        residual_iter.Set(static_cast<InputType>(residual));
                        for (int i = 0; i < ModelType::NV; i++) {
                            output_iters[i].Set(outputs[i]);
                        }
                    }
                    if constexpr (HasDerived) {
                        using DerivedArray = QI_ARRAYN(typename FitType::OutputType, ModelType::ND);
                        DerivedArray derived;
                        m_fit->model.derived(outputs, fixed, derived);
                        for (int i = 0; i < ModelType::ND; i++) {
                            derived_iters[i].Set(derived[i]);
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
                if constexpr (Blocked) {
                    flag_iter.Get().Fill(0);
                    residual_iter.Get().Fill(0);
                    for (auto &o : output_iters) {
                        o.Get().Fill(0);
                    }
                } else {
                    flag_iter.Set(0);
                    flag_iter.Set(0);
                    for (auto &o : output_iters) {
                        o.Set(0);
                    }
                }
                if constexpr (HasDerived) {
                    for (auto &d : derived_iters) {
                        d.Set(0);
                    }
                }
                if (m_allResiduals) {
                    for (auto &r : residuals_iters) {
                        r.Get().Fill(0);
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
            for (int i = 0; i < ModelType::NF; i++) {
                if (this->GetFixed(i))
                    ++fixed_iters[i];
            }
            for (int i = 0; i < ModelType::NV; i++) {
                ++output_iters[i];
            }
            if constexpr (HasDerived) {
                for (auto &d : derived_iters) {
                    ++d;
                }
            }
            ++flag_iter;
            ++residual_iter;
        }
    }
}; // namespace QI

} // namespace QI

#endif // MODELFITFILTER_H
