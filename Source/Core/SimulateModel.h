/*
 *  SimulateModel.h
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QI_SIMULATEMODEL_H
#define QI_SIMULATEMODEL_H

#include "JSON.h"
#include "ModelSimFilter.h"
#include "ImageIO.h"

namespace QI {

template<typename Model, bool MultiOutput>
void SimulateModel(rapidjson::Value &json, const Model &model,
                   const std::vector<std::string> &fixedpaths,const std::vector<std::string> &outpaths,
                   const bool verbose, const double noise) {
    auto simulator = itk::ModelSimFilter<Model, MultiOutput>::New(model);
    simulator->SetNoise(noise);
    QI_LOG(verbose, "Reading varying parameters");
    for (auto i = 0; i < Model::NV; i++) {
        const std::string v = model.varying_names[i];
        const std::string vname = v + "File";
        const std::string vfile = QI::GetMember(json, vname).GetString();
        simulator->SetVarying(i, QI::ReadImage(vfile, verbose));
    }
    if (fixedpaths.size() != Model::NF) {
        QI_FAIL("Number of fixed paths " << fixedpaths.size() << " does not match number of parameters " << Model::NF);
    }
    QI_LOG(verbose, "Reading fixed parameters");
    for (auto i = 0; i < Model::NF; i++) {
        if (fixedpaths[i].size() > 0) {
            simulator->SetFixed(i, QI::ReadImage(fixedpaths[i], verbose));
        }
    }
    QI_LOG(verbose, "Simulating model");
    simulator->Update();
    if constexpr(MultiOutput) {
        if (outpaths.size() != model.NO) {
           QI_FAIL("Number of output paths " << outpaths.size() << " does not match number of outputs " << model.NO);
        }
        for (size_t i = 0; i < model.NO; i++) {
            QI_LOG(verbose, "Writing output image: " << outpaths[i]);
            QI::WriteVectorImage(simulator->GetOutput(i), outpaths[i]);
        }
    } else {
        QI_LOG(verbose, "Writing output image: " << outpaths[0]);
        QI::WriteVectorImage(simulator->GetOutput(0), outpaths[0]);
    }

}

} // End namespace QI

#endif // SIMULATEMODEL_H