/*
 *  transient_main.cpp
 *
 *  Copyright (c) 2019 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#define QI_DEBUG_BUILD 1

#include "Args.h"
#include "Macro.h"
#include "Util.h"

#include "prep_sc.h"

#include <Eigen/Core>

int parmesan_opt(args::Subparser &parser) {
    args::Positional<std::string> json_path(parser, "JSON", "Parameter JSON file");
    // args::Positional<std::string> out_path(parser, "OUTPUT", "Basis JSON file");

    args::ValueFlag<float> tr(parser, "TR", "TR", {"tr"}, 2.e-3);
    args::ValueFlag<float> tramp(parser, "Tramp", "Tramp", {"tramp"}, 10.e-3);
    args::ValueFlag<float> trf(parser, "Trf", "Trf", {"trf"}, 10.e-6);
    args::ValueFlag<int>   sps(parser, "SPS", "SPS", {"sps"}, 64);
    args::ValueFlag<int>   spoils(parser, "SPOILS", "SPOILS", {"spoils"}, 2);
    args::ValueFlag<int>   nP(parser, "NPrep", "Number of preps", {"nprep"}, 8);

    Parse(parser);
    QI::CheckPos(json_path);
    // QI::CheckPos(out_path);
    QI::Info("Reading sequence parameters");

    using ModelType = PrepModel;
    PrepSequence sequence(tr.Get(),
                          tramp.Get(),
                          trf.Get(),
                          0,
                          0,
                          sps.Get(),
                          spoils.Get(),
                          Eigen::ArrayXd::Ones(nP.Get()),
                          Eigen::ArrayXd::Ones(nP.Get()),
                          Eigen::ArrayXd::Ones(nP.Get()) * 1e-3,
                          Eigen::ArrayXd::Zero(nP.Get()));

    ModelType               model{{}, sequence};
    ModelType::VaryingArray v     = (model.lo + model.hi) / 2.;
    v(0)                          = 1; // Set M0 to 1
    ModelType::FixedArray const f = model.fixed_defaults;

    Eigen::MatrixXd FIM(model.NV, model.NV);
    for (int ij = 0; ij < model.NV; ij++) {
        auto const dj = model.dsdθ(v, f, ij);
        for (int ii = 0; ii < model.NV; ii++) {
            auto const di = model.dsdθ(v, f, ii);
            FIM(ii, ij)   = di.matrix().dot(dj.matrix());
        }
    }

    QI_DBMAT(FIM);
    Eigen::MatrixXd const iFIM = FIM.inverse();
    Eigen::VectorXd const CRB  = iFIM.diagonal();
    Eigen::VectorXd const nCRB = CRB.array() * sequence.TR * sequence.SPS / v.square();
    QI_DBMAT(iFIM);
    QI_DBVEC(CRB);
    QI_DBVEC(nCRB);
    // QI_DBVEC(cov);

    QI::Info("Finished.");
    return EXIT_SUCCESS;
}
