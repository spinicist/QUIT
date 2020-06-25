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
void SimulateModel(json &                                    json,
                   Model const &                             model,
                   typename Model::FixedNames const &        fixedpaths,
                   std::array<std::string, Model::NI> const &outpaths,
                   std::string const &                       mask_path,
                   bool const                                verbose,
                   double const                              noise,
                   std::string const &                       subRegion) {
    auto simulator = QI::ModelSimFilter<Model, MultiOutput>::New(model, verbose, subRegion);
    simulator->SetNoise(noise);
    for (auto i = 0; i < Model::NV; i++) {
        const std::string vname = fmt::format("{}_map", model.varying_names[i]);
        const std::string vfile = json.at(vname).get<std::string>();
        QI::Log(verbose, "Reading {} from file: {}", vname, vfile);
        simulator->SetVarying(i, QI::ReadImage(vfile, false));
    }
    if constexpr (Model::NF > 0) {
        if (fixedpaths.size() != Model::NF) {
            QI::Fail("Number of fixed paths {} does not match number of parameters {}",
                     fixedpaths.size(),
                     Model::NF);
        }
        for (auto i = 0; i < Model::NF; i++) {
            if (fixedpaths[i].size() > 0) {
                std::string const &fname = model.fixed_names[i];
                std::string const &ffile = fixedpaths[i];
                QI::Log(verbose, "Reading {} from file: {}", fname, ffile);
                simulator->SetFixed(i, QI::ReadImage(ffile, false));
            }
        }
    }
    if (mask_path != "") {
        simulator->SetMask(QI::ReadImage(mask_path, verbose));
    }
    QI::Log(verbose, "Noise level is {}\nSimulating model...", noise);
    srand((unsigned int)time(0));
    simulator->Update();
    QI::Log(verbose, "Finished");
    if constexpr (MultiOutput) {
        if (outpaths.size() != model.num_outputs()) {
            QI::Fail("Number of output paths {} does not match number of outputs {}",
                     outpaths.size(),
                     model.num_outputs());
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