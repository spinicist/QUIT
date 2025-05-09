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

struct Annealer {
    int                           nK;
    PrepModel::VaryingArray const v;

    auto FIM(PrepSequence const &seq) -> Eigen::MatrixXd {
        PrepModel                   m{{}, seq};
        PrepModel::FixedArray const f = m.fixed_defaults;

        Eigen::MatrixXd FIM(m.NV, m.NV);
        for (int ij = 0; ij < m.NV; ij++) {
            auto const dj = m.dsdθ(v, f, ij);
            for (int ii = 0; ii < m.NV; ii++) {
                auto const di = m.dsdθ(v, f, ii);
                FIM(ii, ij)   = di.matrix().dot(dj.matrix());
            }
        }
        return FIM;
    }

    auto time(PrepSequence const &seq) -> double {
        double const tAcq = (seq.TR * seq.SPS) * seq.preps();
        double const tTotal = tAcq + seq.Tpreseg.sum() + seq.Tpostseg.sum();
        double const time = tTotal / tAcq;
        return time;
    }

    auto nCRB(PrepSequence const &seq) -> Eigen::ArrayXd {
        Eigen::VectorXd const CRB  = FIM(seq).inverse().diagonal();
        Eigen::VectorXd const nCRB = CRB.array() * time(seq) / v.square();
        return nCRB;
    }

    auto cost(PrepSequence const &seq) -> double {
        double cost = nCRB(seq).sum();
        return cost;
    }

    auto perturbation(int N, double T) -> Eigen::ArrayXd {
        return T * (M_PI_2 * Eigen::ArrayXd::Random(N)).tan();
    }

    auto perturb(PrepSequence const &pseq, double const step) -> PrepSequence {
        PrepSequence nseq = pseq;
        int const    N    = pseq.preps();
        nseq.SPS          = std::clamp(int(nseq.SPS + perturbation(1, step * 10)(0)), 32, 512);
        nseq.FA       = (nseq.FA + perturbation(N, step * M_PI / 180)).max(M_PI * 0.5 / 180).min(M_PI * 4. / 180.);
        nseq.FAprep   = (nseq.FAprep + perturbation(N, step * 1 * M_PI / 180)).max(0).min(2 * M_PI);
        nseq.Tpreseg  = (nseq.Tpreseg + perturbation(N, step * 1e-3)).max(0);
        nseq.Tpostseg = (nseq.Tpostseg + perturbation(N, step * 1e-3)).max(0);
        nseq.fprep    = nseq.fprep + perturbation(N, step * 1);
        nseq.Tprep    = (nseq.Tprep + perturbation(N, step * 1e-3)).max(0.5e-3).min(50.e-3);

        nseq.FAprep(0) = M_PI;
        nseq.fprep(0)  = 0;

        return nseq;
    }

    auto run(PrepSequence const &s0) -> PrepSequence {
        PrepSequence s    = s0;
        double       c    = cost(s);
        double const T0   = c / 1e3;
        double       step = 1.;
        for (int ik = 0; ik < nK; ik++) {
            double const T = T0 / (1 + ik);
            for (int ij = 0; ij < 128; ij++) {
                auto const sp = perturb(s, step);
                auto       cp = cost(sp);
                if (cp < c || ((Eigen::ArrayXd::Random(1)(0) + 1.) / 2) < std::exp((c - cp) / T)) {
                    fmt::print(stderr, "k {} T {} New cost {} Time {}\n", ik, T, c, cp, time(sp));
                    c = cp;
                    s = sp;
                    // step /= 2.;
                    break;
                }
            }
        }
        return s;
    }
};

int parmesan_opt(args::Subparser &parser) {
    args::Positional<std::string> json_path(parser, "JSON", "Output JSON file");
    // args::Positional<std::string> out_path(parser, "OUTPUT", "Basis JSON file");

    args::ValueFlag<float> tr(parser, "TR", "TR", {"tr"}, 2.e-3);
    args::ValueFlag<float> tramp(parser, "Tramp", "Tramp", {"tramp"}, 10.e-3);
    args::ValueFlag<float> trf(parser, "Trf", "Trf", {"trf"}, 10.e-6);
    args::ValueFlag<int>   spoils(parser, "SPOILS", "SPOILS", {"spoils"}, 2);
    args::ValueFlag<int>   nP(parser, "NPrep", "Number of preps", {"nprep"}, 8);

    args::ValueFlag<int> k(parser, "k", "Annealing iteration count", {'k', "k"}, 512);

    Parse(parser);
    QI::CheckPos(json_path);
    // QI::CheckPos(out_path);
    QI::Info("Reading sequence parameters");

    using ModelType = PrepModel;
    PrepSequence s0(tr.Get(),
                    tramp.Get(),
                    trf.Get(),
                    128,
                    spoils.Get(),
                    Eigen::ArrayXd::Ones(nP.Get()),
                    Eigen::ArrayXd::Ones(nP.Get()) * 360.,
                    Eigen::ArrayXd::Ones(nP.Get()) * 1e-3,
                    Eigen::ArrayXd::Zero(nP.Get()),
                    Eigen::ArrayXd::Zero(nP.Get()),
                    Eigen::ArrayXd::LinSpaced(nP.Get(), -500., 500.));
    s0.FAprep(0) = M_PI;
    s0.fprep(0)  = 0;

    ModelType               model{{}, s0};
    ModelType::VaryingArray v = (model.lo + model.hi) / 2.;
    v(0)                      = 1; // Set M0 to 1

    Annealer   a{k.Get(), v};
    auto const s = a.run(s0);

    QI_DB(a.nCRB(s));

    auto j = json{{"PrepZTE", s}};
    QI::WriteJSON(json_path.Get(), j);

    QI::Info("Finished.");
    return EXIT_SUCCESS;
}
