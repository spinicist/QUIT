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

#include "ImageIO.h"
#include "JSON.h"
#include "ModelSimFilter.h"

namespace QI {

template <typename Model, bool MultiOutput>
void SimulateModel(rapidjson::Value &json, const Model &model,
                   const std::vector<std::string> &fixedpaths,
                   const std::vector<std::string> &outpaths, const bool verbose,
                   const double noise) {
    auto simulator = QI::ModelSimFilter<Model, MultiOutput>::New(model);
    simulator->SetNoise(noise);
    QI::Log(verbose, "Reading varying parameters");
    for (auto i = 0; i < Model::NV; i++) {
        const std::string v     = model.varying_names[i];
        const std::string vname = v + "File";
        const std::string vfile = QI::GetMember(json, vname).GetString();
        simulator->SetVarying(i, QI::ReadImage(vfile, verbose));
    }
    if (fixedpaths.size() != Model::NF) {
        QI::Fail("Number of fixed paths {} does not match number of parameters {}",
                 fixedpaths.size(), Model::NF);
    }
    QI::Log(verbose, "Reading fixed parameters");
    for (auto i = 0; i < Model::NF; i++) {
        if (fixedpaths[i].size() > 0) {
            simulator->SetFixed(i, QI::ReadImage(fixedpaths[i], verbose));
        }
    }
    QI::Log(verbose, "Simulating model...");
    simulator->Update();
    QI::Log(verbose, "Finished");
    if constexpr (MultiOutput) {
        if (outpaths.size() != model.num_outputs()) {
            QI::Fail("Number of output paths {} does not match number of outputs {}",
                     outpaths.size(), model.num_outputs());
        }
        for (size_t i = 0; i < model.num_outputs(); i++) {
            QI::WriteImage(simulator->GetOutput(i), outpaths[i], verbose);
        }
    } else {
        QI::WriteImage(simulator->GetOutput(0), outpaths[0], verbose);
    }
}

} // End namespace QI

#endif // SIMULATEMODEL_H