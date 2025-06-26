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
#include "ParameterGrid.h"
#include "Util.h"

#include "prep_sc.h"

#include <Eigen/Core>

#include <pagmo/algorithm.hpp>
#include <pagmo/algorithms/cmaes.hpp>
#include <pagmo/algorithms/compass_search.hpp>
#include <pagmo/algorithms/de1220.hpp>
#include <pagmo/algorithms/mbh.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/types.hpp>

namespace {
auto RoundDP(double x, int N) -> double {
    return std::round(x * std::pow(10, N)) / std::pow(10, N);
}

auto B1toFlip(double b1, double t) {
    // ɣ = 42.5777 Hz per uT
    double flip = b1 * 42.577478 * 2 * M_PI * t;
    // fmt::print(stderr, "b1 {} t {} flip {}\n", b1, t, flip);
    return flip;
}

} // namespace

template <typename ModelT> struct PrepProblem {
    PrepSequence    s0;
    bool            gauss;
    double          maxFA, maxB1prep, maxTprep;
    Eigen::ArrayXXd vs;

    using VaryingArray = typename ModelT::VaryingArray;
    using FixedArray   = typename ModelT::FixedArray;

    auto FIM(PrepSequence const &seq, VaryingArray const &v) const -> Eigen::MatrixXd {
        ModelT           m{{}, seq, Eigen::MatrixXd(), gauss};
        FixedArray const f = m.fixed_defaults;

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

    auto nCRB(PrepSequence const &seq, VaryingArray const &v) const -> Eigen::ArrayXd {
        Eigen::VectorXd const CRB  = FIM(seq, v).inverse().diagonal();
        Eigen::VectorXd const nCRB = CRB.array() * time(seq) / v.square();
        return nCRB;
    }

    auto averageCRB(PrepSequence const &s) const -> Eigen::VectorXd {
        VaryingArray c = VaryingArray::Zero();
        for (int ic = 0; ic < vs.cols(); ic++) {
            c += nCRB(s, vs.col(ic));
        }
        c /= vs.size();
        return c;
    }

    auto cost(PrepSequence const &s) const -> double {
        return std::sqrt(averageCRB(s).squaredNorm() / ModelT::NV);
    }

    pagmo::vector_double fitness(const pagmo::vector_double &dv) const {
        PrepSequence s = s0;
        int const    N = s0.preps() - 1;
        for (int ii = 0; ii < N; ii++) {
            s.FA(ii + 1)     = dv[ii] * M_PI / 180.0;
            s.Tprep(ii + 1)  = dv[2 * N + ii];
            s.fprep(ii + 1)  = dv[3 * N + ii];
            s.FAprep(ii + 1) = B1toFlip(dv[N + ii], s.Tprep(ii + 1));
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
        std::fill_n(hi.begin() + N, N, maxB1prep);
        // Tprep
        std::fill_n(lo.begin() + 2 * N, N, 1e-3);
        std::fill_n(hi.begin() + 2 * N, N, maxTprep);
        // fprep
        std::fill_n(lo.begin() + 3 * N, N, -250.0);
        std::fill_n(hi.begin() + 3 * N, N, 250.0);

        return {lo, hi};
    }
};

template <typename ModelT>
void RunCRB(std::vector<std::string> const &jp, int const N, bool const gauss) {
    fmt::print("[{}] [CRB]\n", fmt::join(ModelT::varying_names, ", "));
    // auto const v = QI::RandomPars(ModelT::lo, ModelT::hi, N);
    for (auto const &p : jp) {
        PrepSequence        s0(QI::ReadJSON(p)["PrepZTE"]);
        PrepProblem<ModelT> pp{s0, gauss};
        // pp.vs = v;
        pp.vs = (ModelT::lo + ModelT::hi) / 2;
        fmt::print("[{:4.3E}] {:4.3E} {}\n", fmt::join(pp.averageCRB(s0), ", "), pp.cost(s0), p);
    }
}

int parmesan_crb(args::Subparser &parser) {
    args::PositionalList<std::string> json_paths(parser, "JSON", "Parameter files");
    args::ValueFlag<int>              N(parser, "N", "Number of random samples", {"N", 'N'}, 128);
    args::Flag                        nodf(parser, "F", "No off-resonance", {'f', "nodf"});
    args::Flag                        gauss(parser, "G", "Gauss Prep Pulses", {'g', "gauss"});
    Parse(parser);
    if (nodf) {
        RunCRB<PrepModel2>(json_paths.Get(), N.Get(), gauss);
    } else {
        RunCRB<PrepModel>(json_paths.Get(), N.Get(), gauss);
    }
    QI::Info("Finished.");
    return EXIT_SUCCESS;
}

int parmesan_opt(args::Subparser &parser) {
    args::Positional<std::string> jpath(parser, "JSON", "Output JSON file");
    // args::Positional<std::string> out_path(parser, "OUTPUT", "Basis JSON file");

    args::ValueFlag<float> tr(parser, "TR", "TR", {"tr"}, 2.e-3);
    args::ValueFlag<float> tramp(parser, "Tramp", "Tramp", {"tramp"}, 5.e-3);
    args::ValueFlag<float> trf(parser, "Trf", "Trf", {"trf"}, 10.e-6);
    args::ValueFlag<int>   spoils(parser, "SPOILS", "SPOILS", {"spoils"}, 2);
    args::ValueFlag<int>   nP(parser, "NPrep", "Number of preps", {"nprep"}, 8);

    args::ValueFlag<double> maxFA(parser, "FA", "Maximum FA", {"maxFA"}, 4.);
    args::ValueFlag<double> maxB1prep(parser, "FA", "Maximum B1prep μT", {"maxB1prep"}, 20.);
    args::ValueFlag<double> maxTprep(parser, "T", "Maximum Tprep", {"maxTprep"}, 50.);

    args::ValueFlag<int> N(
        parser, "N", "Number of varying parameters to optimize over", {"N", 'N'}, 128);
    args::ValueFlag<int> nPop(parser, "P", "Population size", {'p', "pop"}, 1024);
    args::ValueFlag<int> nGen(parser, "G", "Number of generations", {'g', "gen"}, 64);

    args::Flag gauss(parser, "G", "Gauss Prep Pulses", {'g', "gauss"});

    Parse(parser);
    QI::CheckPos(jpath);
    // QI::CheckPos(out_path);
    QI::Info("Reading sequence parameters");

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

    PrepProblem<PrepModel> pp{s0, gauss, maxFA.Get(), maxB1prep.Get(), maxTprep.Get()};
    pp.vs = QI::RandomPars(PrepModel::lo, PrepModel::hi, N.Get());
    pagmo::problem    prob{pp};
    pagmo::population pop(prob, nPop.Get());
    pagmo::algorithm  algo = nGen ? pagmo::algorithm{pagmo::cmaes(
                                       nGen.Get(), -1, -1, -1, -1, 0.5, 1e-6, 1e-6, false, true)} :
                                    pagmo::algorithm{pagmo::mbh(pagmo::compass_search(), 16, 0.1)};

    algo.set_verbosity(1);
    pop = algo.evolve(pop);

    int const            Nv = s0.preps() - 1;
    pagmo::vector_double dv = pop.champion_x();
    for (int ii = 0; ii < Nv; ii++) {
        s0.FA(ii + 1)    = RoundDP(dv[ii], 1) * M_PI / 180.0;
        s0.Tprep(ii + 1) = RoundDP(dv[2 * Nv + ii], 3);
        s0.FAprep(ii + 1) =
            RoundDP(B1toFlip(dv[Nv + ii], s0.Tprep(ii + 1)) * 180 / M_PI, 1) * M_PI / 180.0;
        s0.fprep(ii + 1) = RoundDP(dv[3 * Nv + ii], 0);
    }
    auto j = json{{"PrepZTE", s0}};
    QI::WriteJSON(jpath.Get(), j);
    QI::Info("Finished.");
    return EXIT_SUCCESS;
}
