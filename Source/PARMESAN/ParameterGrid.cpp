#include "ParameterGrid.h"

#include "Log.h"

namespace QI {

auto RegularPars(Eigen::ArrayXd const &lo, Eigen::ArrayXd const &hi, Eigen::ArrayXi const &N)
    -> Eigen::ArrayXXd {
    if (lo.size() != hi.size()) {
        QI::Fail("Pars", "Low values had {} elements, hi had {}", lo.size(), hi.size());
    }
    if (hi.size() != N.size()) {
        QI::Fail("Pars", "High values had {} elements, expected {}", hi.size(), N.size());
    }
    auto const nPar = lo.size();

    Eigen::ArrayXd delta(nPar);
    int            nTotal = 1;
    for (int ii = 0; ii < nPar; ii++) {
        if (N[ii] < 1) {
            QI::Fail("Pars", "{} N was less than 1", ii);
        } else if (N[ii] == 1) {
            delta[ii] = 0.f;
        } else {
            delta[ii] = (hi[ii] - lo[ii]) / (N[ii] - 1);
        }
        nTotal *= N[ii];
    }

    Eigen::ArrayXXd p(nPar, nTotal);
    int             ind = 0;

    std::function<void(int, Eigen::ArrayXd)> dimLoop = [&](int dim, Eigen::ArrayXd pars) {
        for (int id = 0; id < N[dim]; id++) {
            pars[dim] = lo[dim] + id * delta[dim];
            if (dim > 0) {
                dimLoop(dim - 1, pars);
            } else {
                p.col(ind++) = pars;
            }
        }
    };
    dimLoop(nPar - 1, Eigen::ArrayXd::Zero(nPar));
    return p;
}

auto RandomPars(Eigen::ArrayXd const &lo, Eigen::ArrayXd const &hi, Eigen::Index const N)
    -> Eigen::ArrayXXd {
    if (lo.size() != hi.size()) {
        QI::Fail("Pars", "Low values had {} elements, hi had {}", lo.size(), hi.size());
    }
    auto const nPar = lo.size();

    Eigen::ArrayXd const range = hi - lo;
    Eigen::ArrayXXd      p(nPar, N);
    for (Eigen::Index ii = 0; ii < N; ii++) {
        p.col(ii) = lo + range * (Eigen::ArrayXd::Random(nPar) * 0.5 + 0.5); // Random is -1 to 1
    }

    return p;
}

} // namespace QI
