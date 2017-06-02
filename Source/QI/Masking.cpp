/*
 *  Masking.cpp
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "QI/Masking.h"

#include <limits>

#include "itkBinaryThresholdImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelShapeKeepNObjectsImageFilter.h"
#include "itkRelabelComponentImageFilter.h"

#include "QI/Macro.h"

namespace QI {

VolumeI::Pointer ThresholdMask(const VolumeF::Pointer &img, const float thresh) {
    typedef itk::BinaryThresholdImageFilter<VolumeF, VolumeI> TThreshFilter;
    auto threshold = TThreshFilter::New();
    threshold->SetInput(img);
    threshold->SetLowerThreshold(thresh);
    threshold->SetUpperThreshold(std::numeric_limits<float>::infinity());
    threshold->SetInsideValue(1);
    threshold->SetOutsideValue(0);
    threshold->Update();
    VolumeI::Pointer mask = threshold->GetOutput();
    mask->DisconnectPipeline();
    return mask;
}

VolumeI::Pointer OtsuMask(const VolumeF::Pointer &img) {
    auto otsuFilter = itk::OtsuThresholdImageFilter<VolumeF, VolumeI>::New();
    otsuFilter->SetInput(img);
    otsuFilter->SetOutsideValue(1);
    otsuFilter->SetInsideValue(0);
    otsuFilter->Update();
    VolumeI::Pointer mask = otsuFilter->GetOutput();
    mask->DisconnectPipeline();
    return mask;
}

VolumeI::Pointer FindLabels(const QI::VolumeI::Pointer &mask, const int size_threshold, int &keep) {
    auto CC = itk::ConnectedComponentImageFilter<QI::VolumeI, QI::VolumeI>::New();
    auto relabel = itk::RelabelComponentImageFilter<QI::VolumeI, QI::VolumeI>::New();
    CC->SetInput(mask);
    relabel->SetInput(CC->GetOutput());
    relabel->Update();
    // Relabel sorts on size by default, so now work out how many make the size threshold
    auto label_sizes = relabel->GetSizeOfObjectsInPixels();
    if (keep > label_sizes.size())
        keep = label_sizes.size();
    for (int i = 0; i < keep; i++) {
        if (label_sizes[i] < size_threshold) {
            keep = i;
            break;
        }
    }
    if (keep == 0) {
        QI_EXCEPTION("No labels found in mask");
    }

    typedef itk::LabelShapeKeepNObjectsImageFilter<QI::VolumeI> TKeepN;
    TKeepN::Pointer keepN = TKeepN::New();
    keepN->SetInput(relabel->GetOutput());
    keepN->SetBackgroundValue(0);
    keepN->SetNumberOfObjects(keep);
    keepN->SetAttribute(TKeepN::LabelObjectType::NUMBER_OF_PIXELS);
    keepN->Update();
    VolumeI::Pointer labels = keepN->GetOutput();
    labels->DisconnectPipeline();
    return labels;
}

} // End namespace QI
