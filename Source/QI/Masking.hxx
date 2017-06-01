/*
 *  Masking.hxx
 *
 *  Copyright (c) 2017 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_MASKING_HXX
#define QI_MASKING_HXX

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

} // End namespace QI

#endif