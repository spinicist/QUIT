#pragma once
#include "ImageTypes.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkMultiThreaderBase.h"

namespace QI {
template <typename ModelType>
void ModelFunc(ModelType const &                                  model,
               std::array<VolumeF::Pointer, ModelType::NI> const &inputs,
               std::array<VolumeF::Pointer, ModelType::NF> const &fixed,
               VolumeF::Pointer const &                           mask,
               itk::ThreadIdType const                            nThreads,
               std::array<VolumeF::Pointer, ModelType::NV> &      outputs) {
    auto mt = itk::MultiThreaderBase::New();
    mt->SetNumberOfWorkUnits(nThreads);
    mt->ParallelizeImageRegion<3>(
        inputs.at(0)->GetBufferedRegion(),
        [&](const QI::VolumeF::RegionType &region) {
            std::array<itk::ImageRegionConstIterator<QI::VolumeF>, ModelType::NI> in_its;
            for (int ii = 0; ii < ModelType::NI; ii++) {
                in_its[ii] = itk::ImageRegionConstIterator<QI::VolumeF>(inputs[ii], region);
            }
            std::array<itk::ImageRegionConstIterator<QI::VolumeF>, ModelType::NF> fixed_its;
            for (int ii = 0; ii < ModelType::NF; ii++) {
                if (fixed[ii]) {
                    fixed_its[ii] = itk::ImageRegionConstIterator<QI::VolumeF>(fixed[ii], region);
                }
            }
            itk::ImageRegionConstIterator<QI::VolumeF> mask_it;
            if (mask) {
                mask_it = itk::ImageRegionConstIterator<QI::VolumeF>(mask, region);
            }
            std::array<itk::ImageRegionIterator<QI::VolumeF>, ModelType::NV> out_its;
            for (int ii = 0; ii < ModelType::NV; ii++) {
                out_its[ii] = itk::ImageRegionIterator<QI::VolumeF>(outputs[ii], region);
            }

            for (in_its[0].GoToBegin(); !in_its[0].IsAtEnd();) {
                if (!mask || mask_it.Get()) {
                    QI_ARRAYN(typename ModelType::DataType, ModelType::NI) in_vals;
                    for (int ii = 0; ii < ModelType::NI; ii++) {
                        in_vals[ii] = in_its[ii].Get();
                    }

                    auto fixed_vals = model.fixed_defaults;
                    for (int ii = 0; ii < ModelType::NF; ii++) {
                        if (fixed[ii]) {
                            fixed_vals[ii] = fixed_its[ii].Get();
                        }
                    }
                    auto const out_vals = model.fit(in_vals, fixed_vals);
                    for (int ii = 0; ii < ModelType::NV; ii++) {
                        out_its[ii].Set(out_vals[ii]);
                    }
                } else {
                    for (int ii = 0; ii < ModelType::NV; ii++) {
                        out_its[ii].Set(0.0);
                    }
                }
                for (int ii = 0; ii < ModelType::NV; ii++) {
                    ++out_its[ii];
                }
                for (int ii = 0; ii < ModelType::NF; ii++) {
                    if (fixed[ii]) {
                        ++fixed_its[ii];
                    }
                }
                for (int ii = 0; ii < ModelType::NI; ii++) {
                    ++in_its[ii];
                }
                if (mask) {
                    ++mask_it;
                }
            }
        },
        nullptr);
}

} // namespace QI