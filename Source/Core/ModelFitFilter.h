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
#include "Model.h"
#include "Monitor.h"
#include "Util.h"

namespace QI {

template <typename FitType>
class ModelFitFilter
    : public itk::ImageToImageFilter<
          itk::VectorImage<typename IOPrecision<typename FitType::ModelType::DataType>::Type, 3>,
          typename BlockTypes<
              FitType::Blocked,
              3,
              typename IOPrecision<typename FitType::ModelType::ParameterType>::Type>::Type> {
  public:
    using ModelType     = typename FitType::ModelType;
    using DataType      = typename ModelType::DataType;
    using ParameterType = typename ModelType::ParameterType;
    using RMSErrorType  = typename FitType::RMSErrorType;

    using VaryingArray  = typename ModelType::VaryingArray;
    using FixedArray    = typename ModelType::FixedArray;
    using CovarArray    = typename ModelType::CovarArray;
    using DataArray     = QI_ARRAY(DataType);
    using ResidualArray = QI_ARRAY(DataType);

    using InputPixelType    = typename IOPrecision<DataType>::Type;
    using OutputPixelType   = typename IOPrecision<ParameterType>::Type;
    using RMSErrorPixelType = typename IOPrecision<RMSErrorType>::Type;
    using FixedPixelType    = typename IOPrecision<ParameterType>::Type;

    static constexpr int ImageDim = 3;

    using TInputImage = itk::VectorImage<InputPixelType, ImageDim>;
    using TFixedImage = itk::Image<FixedPixelType, ImageDim>;
    using TMaskImage  = itk::Image<float, ImageDim>;

    static constexpr bool Blocked = FitType::Blocked;

    using TOutputImage   = typename BlockTypes<Blocked, ImageDim, OutputPixelType>::Type;
    using TFlagImage     = typename BlockTypes<Blocked, ImageDim, typename FitType::FlagType>::Type;
    using TRMSErrorImage = typename BlockTypes<Blocked, ImageDim, RMSErrorPixelType>::Type;
    using TResidualsImage = TInputImage;

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
    static constexpr int DerivedOffset = HasDerived ? ModelType::NV : -1;
    static constexpr int FlagOffset    = HasDerived ? DerivedOffset + ModelType::ND : ModelType::NV;
    static constexpr int RMSErrorOffset  = FlagOffset + 1;
    static constexpr int CovarOffset     = RMSErrorOffset + 1;
    static constexpr int ResidualsOffset = CovarOffset + ModelType::NCov;
    static constexpr int TotalOutputs    = ResidualsOffset + ModelType::NI;

    ModelFitFilter(FitType const *    f,
                   const bool         verbose,
                   const bool         covar,
                   const bool         allResids,
                   std::string const &subregion) :
        m_fit(f),
        m_verbose(verbose), m_allResiduals(allResids), m_covar(covar) {
        this->SetNumberOfRequiredInputs(ModelType::NI);
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
        this->DynamicMultiThreadingOn();
    }

    void SetInput(unsigned int i, const TInputImage *image) override {
        if (static_cast<int>(i) < ModelType::NI) {
            this->SetNthInput(i, const_cast<TInputImage *>(image));
        } else {
            itkExceptionMacro("Requested input " << i << " does not exist (" << ModelType::NI
                                                 << " inputs)");
        }
    }

    typename TInputImage::ConstPointer GetInput(const int i) const {
        if (i < ModelType::NI) {
            return static_cast<const TInputImage *>(this->itk::ProcessObject::GetInput(i));
        } else {
            QI::Fail("Requested input {} but {} has {}", i, typeid(FitType).name(), ModelType::NI);
        }
    }

    void SetFixed(const int i, const TFixedImage *image) {
        if (i < ModelType::NF) {
            this->SetNthInput(ModelType::NI + FixedOffset + i, const_cast<TFixedImage *>(image));
        } else {
            QI::Fail("Tried to set fixed input {} but {} only has {} fixed inputs",
                     i,
                     typeid(ModelType).name(),
                     ModelType::NF);
        }
    }

    typename TFixedImage::ConstPointer GetFixed(const int i) const {
        if (i < ModelType::NF) {
            size_t index = ModelType::NI + FixedOffset + i;
            return static_cast<const TFixedImage *>(this->itk::ProcessObject::GetInput(index));
        } else {
            QI::Fail("Requested fixed input {} but {} has {}",
                     i,
                     typeid(ModelType).name(),
                     ModelType::NF);
        }
    }

    void SetMask(const TMaskImage *mask) {
        this->SetNthInput(ModelType::NI + MaskOffset, const_cast<TMaskImage *>(mask));
    }

    typename TMaskImage::ConstPointer GetMask() const {
        return static_cast<const TMaskImage *>(
            this->itk::ProcessObject::GetInput(ModelType::NI + MaskOffset));
    }

    void SetOutputAllResiduals(const bool r) { m_allResiduals = r; }
    void SetOutputCovar(const bool covar) { m_covar = covar; }

    void SetSubregion(const TRegion &sr) {
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

    TRMSErrorImage *GetRMSErrorOutput() {
        return dynamic_cast<TRMSErrorImage *>(this->itk::ProcessObject::GetOutput(RMSErrorOffset));
    }

    TResidualsImage *GetResidualsOutput(const int i) {
        return dynamic_cast<TResidualsImage *>(
            this->itk::ProcessObject::GetOutput(ResidualsOffset + i));
    }

    TOutputImage *GetCovarOutput(const int i) {
        if (i < ModelType::NCov) {
            return dynamic_cast<TOutputImage *>(
                this->itk::ProcessObject::GetOutput(CovarOffset + i));
        } else {
            QI::Fail("Requested covar output {} but {} has {}",
                     i,
                     typeid(ModelType).name(),
                     ModelType::NV);
        }
    }

    TFlagImage *GetFlagOutput() {
        return dynamic_cast<TFlagImage *>(this->itk::ProcessObject::GetOutput(FlagOffset));
    }

    TOutputImage *GetDerivedOutput(const int i) {
        if constexpr (HasDerived) {
            if (i < ModelType::ND) {
                return dynamic_cast<TOutputImage *>(
                    this->itk::ProcessObject::GetOutput(DerivedOffset + i));
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
    void ReadInputs(std::vector<std::string> const &      inputs,
                    typename ModelType::FixedNames const &fixed,
                    std::string const &                   mask) {
        if (static_cast<size_t>(ModelType::NI) != inputs.size()) {
            QI::Fail("Number of input file paths did not match number of inputs for model");
        }

        for (int i = 0; i < ModelType::NI; i++) {
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
        QI::WriteImage(GetRMSErrorOutput(), prefix + "rmse" + QI::OutExt(), m_verbose);
        QI::WriteImage(GetFlagOutput(), prefix + "iterations" + QI::OutExt(), m_verbose);
        if (m_covar) {
            for (int ii = 0; ii < ModelType::NV; ii++) {
                auto const &name = m_fit->model.varying_names.at(ii);
                QI::WriteImage(
                    GetCovarOutput(ii), prefix + "CoV_" + name + QI::OutExt(), m_verbose);
            }
            int index = ModelType::NV;
            for (int ii = 0; ii < ModelType::NV; ii++) {
                auto const &name1 = m_fit->model.varying_names.at(ii);
                for (int jj = ii + 1; jj < ModelType::NV; jj++) {
                    auto const &name2 = m_fit->model.varying_names.at(jj);
                    QI::WriteImage(GetCovarOutput(index++),
                                   prefix + "Corr_" + name1 + "_" + name2 + QI::OutExt(),
                                   m_verbose);
                }
            }
        }
        if (m_allResiduals) {
            for (int i = 0; i < ModelType::NI; i++) {
                QI::WriteImage(GetResidualsOutput(i),
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
        if (idx < static_cast<itype>(ModelType::NV)) { // Varying parameter outputs
            return TOutputImage::New().GetPointer();
        } else if (idx < static_cast<itype>(FlagOffset)) { // Derived parameter outputs
            return TOutputImage::New().GetPointer();
        } else if (idx == static_cast<itype>(FlagOffset)) {
            return TFlagImage::New().GetPointer();
        } else if (idx == static_cast<itype>(RMSErrorOffset)) {
            return TRMSErrorImage::New().GetPointer();
        } else if (idx < static_cast<itype>(CovarOffset + ModelType::NCov)) {
            return TOutputImage::New().GetPointer();
        } else if (idx < static_cast<itype>(ResidualsOffset + ModelType::NI)) {
            return TResidualsImage::New().GetPointer();
        } else {
            QI::Fail("Attempted to create output {} but {} has {}",
                     idx,
                     typeid(FitType).name(),
                     TotalOutputs);
        }
    }

    const FitType *m_fit;
    const bool     m_verbose, m_allResiduals, m_covar;
    bool           m_hasSubregion = false;
    TRegion        m_subregion;
    int            m_blocks = 1;

    virtual void GenerateOutputInformation() override {
        Superclass::GenerateOutputInformation();

        const auto ip = this->GetInput(0);
        // Verify images are all the same size (ITK checks they have valid orientation)
        for (size_t i = 1; i < this->GetNumberOfRequiredInputs(); i++) {
            const auto ip2 = this->GetInput(i);
            if (ip->GetLargestPossibleRegion() != ip2->GetLargestPossibleRegion()) {
                QI::Fail("Input parameter images are not all the same size");
            }
        }

        for (int i = 0; i < ModelType::NI; i++) {
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

        auto f = this->GetFlagOutput();
        f->SetRegions(region);
        f->SetSpacing(spacing);
        f->SetOrigin(origin);
        f->SetDirection(direction);
        if constexpr (Blocked) {
            f->SetNumberOfComponentsPerPixel(m_blocks);
        }
        f->Allocate(true);

        auto rms = this->GetRMSErrorOutput();
        rms->SetRegions(region);
        rms->SetSpacing(spacing);
        rms->SetOrigin(origin);
        rms->SetDirection(direction);
        if constexpr (Blocked) {
            rms->SetNumberOfComponentsPerPixel(m_blocks);
        }
        rms->Allocate(true);

        if (m_covar) {
            for (int ii = 0; ii < ModelType::NCov; ii++) {
                auto op = this->GetCovarOutput(ii);
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
            for (int i = 0; i < ModelType::NI; i++) {
                auto res = this->GetResidualsOutput(i);
                res->SetRegions(region);
                res->SetSpacing(spacing);
                res->SetOrigin(origin);
                res->SetDirection(direction);
                res->SetNumberOfComponentsPerPixel(m_fit->input_size(i) * m_blocks);
                res->Allocate(true);
            }
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

        std::vector<itk::ImageRegionConstIterator<TInputImage>> input_iters(ModelType::NI);
        std::vector<itk::ImageRegionIterator<TResidualsImage>>  residuals_iters(ModelType::NI);
        for (int i = 0; i < ModelType::NI; i++) {
            input_iters[i] = itk::ImageRegionConstIterator<TInputImage>(this->GetInput(i), region);
            if (m_allResiduals) {
                residuals_iters[i] =
                    itk::ImageRegionIterator<TResidualsImage>(this->GetResidualsOutput(i), region);
            }
        }
        std::array<itk::ImageRegionConstIterator<TFixedImage>, ModelType::NF> fixed_iters;
        for (int i = 0; i < ModelType::NF; i++) {
            typename TFixedImage::ConstPointer c = this->GetFixed(i);
            if (c) {
                fixed_iters[i] = itk::ImageRegionConstIterator<TFixedImage>(c, region);
            }
        }

        std::array<itk::ImageRegionIterator<TOutputImage>, ModelType::NV> output_iters;
        for (int i = 0; i < ModelType::NV; i++) {
            output_iters[i] = itk::ImageRegionIterator<TOutputImage>(this->GetOutput(i), region);
        }

        std::array<itk::ImageRegionIterator<TOutputImage>, ModelType::ND> derived_iters;
        if constexpr (HasDerived) {
            for (int i = 0; i < ModelType::ND; i++) {
                derived_iters[i] =
                    itk::ImageRegionIterator<TOutputImage>(this->GetDerivedOutput(i), region);
            }
        }

        std::array<itk::ImageRegionIterator<TOutputImage>, ModelType::NCov> covar_iters;
        if (m_covar) {
            for (int ii = 0; ii < ModelType::NCov; ii++) {
                covar_iters[ii] =
                    itk::ImageRegionIterator<TOutputImage>(this->GetCovarOutput(ii), region);
            }
        }

        // Keep the index on this one to report voxel locations where algorithm fails.
        itk::ImageRegionIteratorWithIndex<TRMSErrorImage> rmse_iter(this->GetRMSErrorOutput(),
                                                                    region);
        itk::ImageRegionIterator<TFlagImage>              flag_iter(this->GetFlagOutput(), region);

        VaryingArray outputs;
        FixedArray   fixed;
        CovarArray * covar = m_covar ? new CovarArray : nullptr;

        while (!input_iters[0].IsAtEnd()) {
            if (!mask || mask_iter.Get()) {
                for (int b = 0; b < m_blocks; b++) {
                    std::vector<DataArray> inputs(ModelType::NI);
                    for (int i = 0; i < ModelType::NI; i++) {
                        auto input_data       = input_iters[i].Get();
                        inputs[i]             = DataArray(m_fit->input_size(i));
                        const int block_start = b * m_fit->input_size(i);
                        for (Eigen::Index j = 0; j < m_fit->input_size(i); j++) {
                            inputs[i][j] = input_data[j + block_start];
                        }
                    }

                    outputs = VaryingArray::Zero();
                    if constexpr (ModelType::NF > 0) {
                        fixed = m_fit->model.fixed_defaults;
                        for (size_t i = 0; i < fixed_iters.size(); i++) {
                            if (this->GetFixed(i)) {
                                fixed[i] = fixed_iters[i].Get();
                            }
                        }
                    }
                    if (m_covar) {
                        *covar = CovarArray::Zero();
                    }

                    typename FitType::RMSErrorType rmse = 0;
                    typename FitType::FlagType     flag = 0;
                    std::vector<ResidualArray>     rs; // Leave size 0 if user doesn't want them
                    if (m_allResiduals) {
                        for (int i = 0; i < ModelType::NI; i++) {
                            const int block_size = m_fit->input_size(i);
                            rs.push_back(ResidualArray::Zero(block_size));
                        }
                    }

                    QI::FitReturnType status;
                    if constexpr (Blocked && Indexed) {
                        status = m_fit->fit(
                            inputs, fixed, outputs, covar, rmse, rs, flag, b, rmse_iter.GetIndex());
                    } else if constexpr (Blocked) {
                        status = m_fit->fit(inputs, fixed, outputs, covar, rmse, rs, flag, b);
                    } else if constexpr (Indexed) {
                        status = m_fit->fit(
                            inputs, fixed, outputs, covar, rmse, rs, flag, rmse_iter.GetIndex());
                    } else {
                        status = m_fit->fit(inputs, fixed, outputs, covar, rmse, rs, flag);
                    }

                    if (!status.success && m_verbose) {
                        QI::Warn(
                            "Fit failed for voxel {}: {}", rmse_iter.GetIndex(), status.message);
                    }

                    if constexpr (Blocked) {
                        flag_iter.Get()[b] = flag;
                        rmse_iter.Get()[b] = rmse;
                        for (int i = 0; i < ModelType::NV; i++) {
                            output_iters[i].Get()[b] = outputs[i];
                        }
                    } else {
                        flag_iter.Set(flag);
                        rmse_iter.Set(rmse);
                        for (int i = 0; i < ModelType::NV; i++) {
                            output_iters[i].Set(outputs[i]);
                        }
                        if (m_covar) {
                            for (int ii = 0; ii < ModelType::NCov; ii++) {
                                covar_iters[ii].Set((*covar)[ii]);
                            }
                        }
                    }
                    if constexpr (HasDerived) {
                        typename ModelType::DerivedArray derived;
                        m_fit->model.derived(outputs, fixed, derived);
                        for (int i = 0; i < ModelType::ND; i++) {
                            derived_iters[i].Set(derived[i]);
                        }
                    }
                    if (m_allResiduals) {
                        for (int i = 0; i < ModelType::NI; i++) {
                            const int block_start = m_fit->input_size(i) * b;
                            for (int j = 0; j < m_fit->input_size(i); j++) {
                                residuals_iters[i].Get()[j + block_start] = rs[i][j];
                            }
                        }
                    }
                }
            } else {
                if constexpr (Blocked) {
                    flag_iter.Get().Fill(0);
                    rmse_iter.Get().Fill(0);
                    for (auto &o : output_iters) {
                        o.Get().Fill(0);
                    }
                } else {
                    flag_iter.Set(0);
                    flag_iter.Set(0);
                    for (auto &o : output_iters) {
                        o.Set(0);
                    }
                    if (m_covar) {
                        for (auto &c : covar_iters) {
                            c.Set(0);
                        }
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
            for (int i = 0; i < ModelType::NI; i++) {
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
            if (m_covar) {
                for (auto &c : covar_iters) {
                    ++c;
                }
            }
            ++flag_iter;
            ++rmse_iter;
        }
    }
}; // namespace QI

} // namespace QI

#endif // MODELFITFILTER_H
