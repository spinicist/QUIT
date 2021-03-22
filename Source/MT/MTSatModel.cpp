#include "MTSatModel.h"

auto MTSatModel::fit(QI_ARRAYN(DataType, NI) const &     in,
                     QI_ARRAYN(ParameterType, NF) const &fixed) const
    -> QI_ARRAYN(ParameterType, NV) {
    auto const  S_pd = in[0];
    auto const  S_t1 = in[1];
    auto const  S_mt = in[2];
    auto        B1   = fixed[0];
    auto const &s    = sequence;

    auto const R1 =
        std::clamp((B1 * B1 / 2.) * (S_t1 * s.al_t1 / s.TR_t1 - S_pd * s.al_pd / s.TR_pd) /
                       (S_pd / s.al_pd - S_t1 / s.al_t1),
                   0.,
                   10.);
    auto const A =
        std::max((S_pd * S_t1 / B1) * (s.TR_pd * s.al_t1 / s.al_pd - s.TR_t1 * s.al_pd / s.al_t1) /
                     (S_t1 * s.TR_pd * s.al_t1 - S_pd * s.TR_t1 * s.al_pd),
                 0.);
    auto const d           = (A * s.al_mt / S_mt - 1.0) * R1 * s.TR_mt - s.al_mt * s.al_mt / 2;
    auto const d_corrected = std::clamp(d * (1.0 - C) / (1.0 - C * B1), 0., 0.1);
    return {A, R1, d_corrected * 100};
}