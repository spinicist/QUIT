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

// #define QI_DEBUG_BUILD 1

#include "Args.h"
#include "Macro.h"
#include "Util.h"

#include "prep_sc.h"

#include <Eigen/Core>

#include <pagmo/algorithm.hpp>
#include <pagmo/algorithms/de1220.hpp>
#include <pagmo/algorithms/mbh.hpp>
#include <pagmo/algorithms/nsga2.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/types.hpp>

struct PrepProblem {
    PrepSequence                         s0;
    double                               maxFA, maxFAprep, maxTprep;
    std::vector<PrepModel::VaryingArray> vs;

    auto FIM(PrepSequence const &seq, PrepModel::VaryingArray const &v) const -> Eigen::MatrixXd {
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

    auto time(PrepSequence const &seq) const -> double {
        double const tAcq   = (seq.TR * seq.SPS + 2 * seq.Tramp) * seq.preps();
        double const tTotal = tAcq + seq.Tprep.sum() + seq.Tpreseg.sum() + seq.Tpostseg.sum();
        return tTotal;
    }

    auto nCRB(PrepSequence const &seq, PrepModel::VaryingArray const &v) const -> Eigen::ArrayXd {
        Eigen::VectorXd const CRB  = FIM(seq, v).inverse().diagonal();
        Eigen::VectorXd const nCRB = CRB.array() * time(seq) / v.square();
        return nCRB;
    }

    auto cost(PrepSequence const &s) const -> double {
        double sum = 0;
        for (auto const &v : vs) {
            double const c = nCRB(s, v).square().sum();
            sum += c;
        }
        return std::sqrt(sum);
    }

    pagmo::vector_double fitness(const pagmo::vector_double &dv) const {
        PrepSequence s = s0;
        int const    N = s0.preps() - 1;
        for (int ii = 0; ii < N; ii++) {
            s.FA(ii + 1)     = dv[ii] * M_PI / 180.0;
            s.FAprep(ii + 1) = dv[N + ii] * M_PI / 180.0;
            s.Tprep(ii + 1)  = dv[2 * N + ii];
            s.fprep(ii + 1)  = dv[3 * N + ii];
        }

        return {cost(s)};
    }

    auto get_nec() const -> pagmo::vector_double::size_type { return 0; }
    auto get_nic() const -> pagmo::vector_double::size_type { return 0; }

    std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const {
        int const N  = s0.preps() - 1;
        auto      lo = pagmo::vector_double(N * 4);
        auto      hi = pagmo::vector_double(N * 4);

        // FA
        std::fill_n(lo.begin(), N, 0.0);
        std::fill_n(hi.begin(), N, maxFA);
        // FAprep
        std::fill_n(lo.begin() + N, N, 0.0);
        std::fill_n(hi.begin() + N, N, maxFAprep);
        // Tprep
        std::fill_n(lo.begin() + 2 * N, N, 1e-3);
        std::fill_n(hi.begin() + 2 * N, N, maxTprep);
        // fprep
        std::fill_n(lo.begin() + 3 * N, N, -250.0);
        std::fill_n(hi.begin() + 3 * N, N, 250.0);

        return {lo, hi};
    }
};

int parmesan_crb(args::Subparser &parser) {
    args::Positional<std::string> json_path(parser, "JSON", "Parameter file");
    Parse(parser);
    QI::CheckPos(json_path);
    PrepSequence            s0(QI::ReadJSON(json_path.Get())["PrepZTE"]);
    PrepProblem             pp{s0};
    PrepModel::VaryingArray v;
    v << 1.0, 1.0, 0.1, 1.0, 10.0; // M0, T1, T2, B1, f0
    pp.vs.push_back(v);
    v << 1.0, 1.0, 0.1, 1.0, -50.0;
    pp.vs.push_back(v);
    v << 1.0, 1.0, 0.1, 1.0, 50.0;
    pp.vs.push_back(v);

    fmt::print("CRB {:4.3f} Time {:4.3f}\n", pp.cost(s0), pp.time(s0));
    fmt::print("[{}]\n", fmt::join(PrepModel::varying_names, ", "));
    fmt::print("[{:4.3f}]\n", fmt::join(pp.nCRB(s0, pp.vs.front()), ", "));
    

    QI::Info("Finished.");
    return EXIT_SUCCESS;
}

int parmesan_opt(args::Subparser &parser) {
    args::Positional<std::string> json_path(parser, "JSON", "Output JSON file");
    // args::Positional<std::string> out_path(parser, "OUTPUT", "Basis JSON file");

    args::ValueFlag<float> tr(parser, "TR", "TR", {"tr"}, 2.e-3);
    args::ValueFlag<float> tramp(parser, "Tramp", "Tramp", {"tramp"}, 5.e-3);
    args::ValueFlag<float> trf(parser, "Trf", "Trf", {"trf"}, 10.e-6);
    args::ValueFlag<int>   spoils(parser, "SPOILS", "SPOILS", {"spoils"}, 2);
    args::ValueFlag<int>   nP(parser, "NPrep", "Number of preps", {"nprep"}, 8);

    args::ValueFlag<double> maxFA(parser, "FA", "Maximum FA", {"maxFA"}, 4.);
    args::ValueFlag<double> maxFAprep(parser, "FA", "Maximum FA", {"maxFAprep"}, 90.);
    args::ValueFlag<double> maxTprep(parser, "T", "Maximum Tprep", {"maxTprep"}, 50.);

    args::ValueFlag<int> nPop(parser, "P", "Population size", {'p', "pop"}, 1024);
    args::ValueFlag<int> nGen(parser, "G", "Number of generations", {'g', "gen"}, 64);

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
                    Eigen::ArrayXd::Ones(nP.Get()) * 4,
                    Eigen::ArrayXd::Ones(nP.Get()) * 80.,
                    Eigen::ArrayXd::Ones(nP.Get()) * 20e-3,
                    Eigen::ArrayXd::Zero(nP.Get()),
                    Eigen::ArrayXd::Zero(nP.Get()),
                    Eigen::ArrayXd::LinSpaced(nP.Get(), -100, 100.));
    s0.FA(0)     = 4 * M_PI / 180.;
    s0.FAprep(0) = M_PI;
    s0.fprep(0)  = 0;
    s0.Tprep(0)  = 1e-3;

    PrepProblem             pp{s0, maxFA.Get(), maxFAprep.Get(), maxTprep.Get()};
    ModelType::VaryingArray v;
    v << 1.0, 1.0, 0.1, 1.0, 10.0; // M0, T1, T2, B1, f0
    pp.vs.push_back(v);
    v << 1.0, 1.0, 0.1, 1.0, -50.0;
    pp.vs.push_back(v);
    v << 1.0, 1.0, 0.1, 1.0, 50.0;
    pp.vs.push_back(v);
    pagmo::problem    prob{pp};
    pagmo::population pop(prob, nPop.Get());
    pagmo::algorithm  algo{pagmo::mbh()};

    algo.set_verbosity(1);
    pop = algo.evolve(pop);

    int const            N  = s0.preps() - 1;
    pagmo::vector_double dv = pop.champion_x();
    for (int ii = 0; ii < N; ii++) {
        s0.FA(ii + 1)     = dv[ii] * M_PI / 180.0;
        s0.FAprep(ii + 1) = dv[N + ii] * M_PI / 180.0;
        s0.Tprep(ii + 1)  = dv[2 * N + ii];
        s0.fprep(ii + 1)  = dv[3 * N + ii];
    }
    auto j = json{{"PrepZTE", s0}};
    QI::WriteJSON(json_path.Get(), j);
    QI::Info("Finished.");
    return EXIT_SUCCESS;
}
