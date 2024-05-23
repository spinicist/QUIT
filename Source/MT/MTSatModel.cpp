#include "MTSatModel.h"

auto MTSatModel::fit(QI_ARRAYN(DataType, NI) const &     in,
                     QI_ARRAYN(ParameterType, NF) const &fixed) const
    -> QI_ARRAYN(ParameterType, NV) {
    auto const  S_pd = in[0];
    auto const  S_t1 = in[1];
    auto const  S_mt = in[2];
    auto        B1   = fixed[0];
    auto const &s    = sequence;
    
    auto al_t1_b1corr = B1 * s.al_t1;
    auto al_pd_b1corr = B1 * s.al_pd;
    if (!smallangle) {
        al_t1_b1corr = 2 * tan(al_t1_b1corr / 2);
        al_pd_b1corr = 2 * tan(al_pd_b1corr / 2);
    }
    
    auto const R1 =
        std::clamp(0.5 * (S_t1 * al_t1_b1corr / s.TR_t1 - S_pd * al_pd_b1corr / s.TR_pd) /
                       (S_pd / al_pd_b1corr - S_t1 / al_t1_b1corr),
                   0.,
                   r1_max);
    auto const A =
        std::max((S_pd * S_t1) * (s.TR_pd * al_t1_b1corr / al_pd_b1corr - s.TR_t1 * al_pd_b1corr / al_t1_b1corr) /
                     (S_t1 * s.TR_pd * al_t1_b1corr - S_pd * s.TR_t1 * al_pd_b1corr),
                 0.);
    auto const d           = (A * s.al_mt / S_mt - 1.0) * R1 * s.TR_mt - s.al_mt * s.al_mt / 2;
    auto const d_corrected = std::clamp(d * (1.0 - C) / (1.0 - C * B1), 0., delta_max / 100.);
    return {A, R1, d_corrected * 100};
}